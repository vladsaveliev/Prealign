from collections import OrderedDict, defaultdict
from itertools import dropwhile
import re
import os
from os.path import join, isfile, isdir, basename, exists, dirname, realpath
import shutil
import traceback

from ngs_utils.logger import critical, err, info, warn, debug
from ngs_utils.file_utils import verify_dir, verify_file, splitext_plus, safe_mkdir, file_transaction, can_reuse


def _sample_name_special_chars(sn):
    fixed_sn = re.sub(r'[\W_]+', r'[\W_]+', sn)
    if fixed_sn.endswith('+'):
        fixed_sn = fixed_sn[:-1] + '*'
    return fixed_sn

def get_hiseq4000_miseq_regexp(sample, suf):
    # sn = ''.join(c for c in sample.name if c.isalnum() or c in ['-', '_', '.', '/'])
    return _sample_name_special_chars(sample.name) + '_S\d+_L\d\d\d_' + suf + '.*\.fastq\.gz'

def get_hiseq_regexp(sample, suf):
    # sn = ''.join(c for c in sample.name if c.isalnum() or c in ['-', '_', '.', '/'])
    return _sample_name_special_chars(sample.name) + '_' + sample.index + '_L\d\d\d_' + suf + '.*\.fastq\.gz'

def get_nextseq500_regexp(sample, suf):
    # sn = ''.join(c for c in sample.name if c.isalnum() or c in ['-', '_', '.', '/'])
    return _sample_name_special_chars(sample.name) + '_S\d+_' + suf + '.*\.fastq\.gz'


class DatasetStructure:
    pre_fastqc_repr = 'Preproc FastQC'
    downsample_targqc_repr = 'TargQC downsampled'

    @staticmethod
    def create(input_dir, proj_infos, samplesheet=None, **kwargs):
        if 'datasets/miseq/' in input_dir.lower():
            return MiSeqStructure(input_dir, proj_infos, samplesheet, **kwargs)

        elif 'datasets/hiseq/' in input_dir.lower():
            return HiSeqStructure(input_dir, proj_infos, samplesheet, **kwargs)

        elif 'datasets/hiseq4000/' in input_dir.lower():
            return HiSeq4000Structure(input_dir, proj_infos, samplesheet, **kwargs)

        elif 'datasets/nextseq500' in input_dir.lower():
            return NextSeq500Structure(input_dir, proj_infos, samplesheet, **kwargs)

        else:
            critical('Directory must be datasets/miseq/, datasets/nextseq500, datasets/hiseq/, or datasets/hiseq4000/. Found ' + input_dir)

    def __init__(self, input_dir, proj_infos, samplesheet=None, **kwargs):
        self.proj_infos = proj_infos

        illumina_project_name = None
        if '/Unalign' in input_dir:
            self.illumina_dir = input_dir.split('/Unalign')[0]
            self.unaligned_dirpath = self.__find_unaligned_dir()
            verify_dir(self.unaligned_dirpath, description='Unalign dir', is_critical=True)
            illumina_project_name = input_dir.split('/Unalign')[1].split('/', 1)[1]  # something like AURA.FFPE.AZ300, in contast with project_name which is something like Bio_123_AURA_FFPE_AZ300
            info('Processing sub-project ' + illumina_project_name)
        else:
            self.illumina_dir = input_dir
            self.unaligned_dirpath = self.__find_unaligned_dir()

        self.basecalls_dirpath = join(self.illumina_dir, 'Data/Intensities/BaseCalls')
        verify_dir(self.basecalls_dirpath, is_critical=True)

        self.basecalls_reports_dirpath = None
        self.bcl2fastq_dirpath = None
        self.source_fastq_dirpath = None

        if samplesheet:
            self.samplesheet_fpath = samplesheet
        else:
            self.samplesheet_fpath = self.__find_sample_sheet()
        self.project_by_name = self._parse_sample_sheet(self.samplesheet_fpath)

        if illumina_project_name:  # we want only a specific project
            if illumina_project_name not in self.project_by_name:
                info()
                critical('Err: project ' + illumina_project_name + ' not in the SampleSheet ' + self.samplesheet_fpath)
            else:
                self.project_by_name = {illumina_project_name: self.project_by_name[illumina_project_name]}

    def __find_unaligned_dir(self):
        unaligned_dirpath = join(self.illumina_dir, 'Unalign')
        if verify_dir(unaligned_dirpath, description='"Unalign" directory', silent=True):
            unaligned_dirpath = unaligned_dirpath
        else:
            unaligned_dirpath = None
            warn('No unalign directory')
        return unaligned_dirpath

    def __find_sample_sheet(self):
        ss_fpath = join(self.illumina_dir, 'SampleSheet.csv')
        if not isfile(ss_fpath):
            ss_fpath = join(self.basecalls_dirpath, 'SampleSheet.csv')
        verify_file(ss_fpath, description='Sample sheet', is_critical=True)
        return ss_fpath

    def _parse_sample_sheet(self, sample_sheet_fpath):
        info('Parsing sample sheet ' + sample_sheet_fpath)
        with open(sample_sheet_fpath) as f:
            def check_if_header(l):
                return any(l.startswith(w) for w in [
                    'Sample_ID,',  # MiSeq
                    'FCID,',       # HiSeq
                    'Lane,',       # HiSeq4000
                ])

            sample_lines = dropwhile(lambda l: not check_if_header(l), f)
            sample_infos = []
            keys = []
            for l in sample_lines:
                if check_if_header(l):
                    keys = l.strip().split(',')
                else:
                    fs = l.strip().split(',')
                    sample_infos.append(dict(zip(keys, fs)))

        project_by_name = OrderedDict()

        for i, info_d in enumerate(sample_infos):
            proj_name = info_d.get('Sample_Project', info_d.get('SampleProject', info_d.get('Project')))
            if proj_name is None:
                warn('  no SampleProject or Sample_Project or Project field in the SampleSheet ' + sample_sheet_fpath)
            elif not proj_name:
                warn('  SampleProject/Sample_Project/Project field is empty in the SampleSheet ' + sample_sheet_fpath)
                     # ', using ' + self.az_prjname_by_subprj[''])
                # proj_name = self.az_prjname_by_subprj['']
            proj_name = proj_name or ''
            if proj_name is not None and proj_name not in project_by_name:
                project_by_name[proj_name] = DatasetProject(proj_name)
            project = project_by_name[proj_name]

            sname = info_d.get('Sample_Name') or info_d.get('SampleName') or info_d.get('SampleRef')
            if sname in project.sample_by_name:
                s = project.sample_by_name[sname]
                s.lane_numbers.add(info_d.get('Lane', 1))  # lanes are in HiSeq and HiSeq4000 (not in MiSeq!)
            else:
                s = DatasetSample(sname, index=info_d.get('index', info_d.get('Index')))
                info('  ' + proj_name + ': ' + s.name)
                s.lane_numbers.add(info_d.get('Lane', 1))  # lanes are in HiSeq and HiSeq4000 (not in MiSeq!)
                if 'FCID' in info_d:
                    s.fcid = info_d['FCID']  # HiSeq only
                project.sample_by_name[sname] = s
                # sample_names.append(info_d[key].replace(' ', '-') + '_' + info_d['Index'] + '_L%03d' % lane)
                # sample_names.append(info_d[key].replace(' ', '-').replace('_', '-') + '_S' + str(i + 1) + '_L001')

        return project_by_name

    def _get_output_dir(self, proj_infos, pname):
        if len(proj_infos) != len(self.project_by_name):
            critical('The number of projects in the input config != the number of projects SampleSheet')

        analysis_dir = None
        output_dir = None
        az_proj_name = None
        if len(proj_infos) == 1:
            analysis_dir = list(proj_infos.values())[0].analysis_dir
            output_dir = list(proj_infos.values())[0].output_dir
            az_proj_name = list(proj_infos.values())[0].project_name
        elif pname in proj_infos:
            analysis_dir = proj_infos[pname].analysis_dir
            output_dir = proj_infos[pname].output_dir
            az_proj_name = proj_infos[pname].project_name
        else:
            critical('Error: cannot correspond the subproject name in the SampleSheet (' + pname + ') and the lines in the conf. ' +
                 'Please, follow the SOP for multiple-project run: http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Pre+Processing+QC+Reporting')
        safe_mkdir(output_dir)
        return analysis_dir, output_dir, az_proj_name


class HiSeqStructure(DatasetStructure):
    def __init__(self, input_dir, proj_infos, samplesheet=None, **kwargs):
        info('Parsing the HiSeq project structure')
        self.kind = 'hiseq'
        DatasetStructure.__init__(self, input_dir, proj_infos, samplesheet=samplesheet, **kwargs)

        verify_dir(self.unaligned_dirpath, is_critical=True)

        self.basecall_stat_html_reports = self.__get_basecall_stats_reports()

        for pname, project in self.project_by_name.items():
            ds_proj_dir = join(self.unaligned_dirpath, 'Project_' + pname.replace(' ', '-'))  #.replace('-', '_').replace('.', '_'))

            analysis_dir, output_dir, az_proj_name = self._get_output_dir(proj_infos, pname)

            project.set_dirpath(ds_proj_dir, analysis_dir, output_dir, az_proj_name)
            for sname, sample in project.sample_by_name.items():
                sample.source_fastq_dirpath = join(project.ds_dir, 'Sample_' + sname.replace(' ', '-'))  #.replace('-', '_').replace('.', '_'))
                sample.set_up_out_dirs(project.fastq_dirpath, project.fastqc_dirpath, project.downsample_targqc_dirpath)

            basecalls_symlink = join(project.ds_dir, 'BaseCallsReports')
            if not exists(basecalls_symlink):
                info('Creating BaseCalls symlink ' + self.basecalls_dirpath + ' -> ' + basecalls_symlink)
                try:
                    os.symlink(self.basecalls_dirpath, basecalls_symlink)
                except OSError:
                    err('Cannot create symlink')
                    traceback.print_exc()
                else:
                    info('Created')
            if exists(basecalls_symlink):
                self.basecalls_dirpath = basecalls_symlink

        self.get_fastq_regexp_fn = get_hiseq_regexp

    # def __get_bcl2fastq_dirpath(self):
    #     # Reading project name
    #     bcl2fastq_dirpath = None
    #     try:
    #         bcl2fastq_dirpath = join(self.unaligned_dirpath, next(fn for fn in os.listdir(self.unaligned_dirpath) if fn.startswith('Project_')))
    #     except StopIteration:
    #         critical('Could not find directory starting with Project_ in ' + self.unaligned_dirpath)
    #     return bcl2fastq_dirpath

    def __get_basecall_stats_reports(self):
        basecall_stats_dirnames = [fname for fname in os.listdir(self.unaligned_dirpath) if fname.startswith('Basecall_Stats_')]
        if len(basecall_stats_dirnames) > 1:
            err('More than 1 Basecall_Stats_* dirs found in unalign_dirpath')
        if len(basecall_stats_dirnames) == 0:
            err('No Basecall_Stats_* dirs found in unalign_dirpath')
        if len(basecall_stats_dirnames) == 1:
            self.basecalls_reports_dirpath = join(self.unaligned_dirpath, basecall_stats_dirnames[0])
            basecall_reports = [verify_file(join(self.basecalls_reports_dirpath, html_fname)) for html_fname in
                                ['Demultiplex_Stats.htm', 'All.htm', 'IVC.htm']]
            return filter(None, basecall_reports)


class MiSeqStructure(DatasetStructure):
    def __init__(self, input_dir, proj_infos, samplesheet=None, **kwargs):
        info('Parsing the MiSeq project structure')
        self.kind = 'miseq'
        DatasetStructure.__init__(self, input_dir, proj_infos, samplesheet=samplesheet, **kwargs)

        base_dirpath = self.unaligned_dirpath
        if not verify_dir(base_dirpath, silent=True):
            base_dirpath = self.basecalls_dirpath
        verify_dir(base_dirpath, description='Source fastq dir')

        for pname, project in self.project_by_name.items():
            ds_proj_dir = join(base_dirpath, pname)
            if not verify_dir(ds_proj_dir, silent=True):
                ds_proj_dir = base_dirpath

            analysis_dir, output_dir, az_proj_name = self._get_output_dir(proj_infos, pname)

            project.set_dirpath(ds_proj_dir, analysis_dir, output_dir, az_proj_name)
            for sample in project.sample_by_name.values():
                sample.source_fastq_dirpath = project.ds_dir
                sample.set_up_out_dirs(project.fastq_dirpath, project.fastqc_dirpath, project.downsample_targqc_dirpath)

        self.basecall_stat_html_reports = self.__get_basecall_stats_reports()
        if self.basecall_stat_html_reports:
            info('basecall_stat_html_reports: ' + str(self.basecall_stat_html_reports))

        self.get_fastq_regexp_fn = get_hiseq4000_miseq_regexp

    def __find_fastq_dir(self):
        for dname in os.listdir(self.unaligned_dirpath):
            dpath = join(self.unaligned_dirpath, dname)
            if isdir(dpath) and any(f.endswith('.fastq.gz') for f in os.listdir(dpath)):
                return dpath

    def __get_basecall_stats_reports(self):
        self.basecalls_reports_dirpath = join(self.unaligned_dirpath, 'Reports', 'html')
        index_html_fpath = join(self.basecalls_reports_dirpath, 'index.html')
        if verify_dir(self.illumina_dir) and verify_file(index_html_fpath):
            return [index_html_fpath]


class HiSeq4000Structure(DatasetStructure):
    def __init__(self, input_dir, proj_infos, samplesheet=None, **kwargs):
        info('Parsing the HiSeq4000 project structure')
        self.kind = 'hiseq4000'
        DatasetStructure.__init__(self, input_dir, proj_infos, samplesheet=samplesheet, **kwargs)

        verify_dir(self.unaligned_dirpath, is_critical=True)

        for pname, project in self.project_by_name.items():
            ds_proj_dir = join(self.unaligned_dirpath, pname)

            analysis_dir, output_dir, az_proj_name = self._get_output_dir(proj_infos, pname)

            project.set_dirpath(ds_proj_dir, analysis_dir, output_dir, az_proj_name)
            for sample in project.sample_by_name.values():
                sample.source_fastq_dirpath = project.ds_dir
                sample.set_up_out_dirs(project.fastq_dirpath, project.fastqc_dirpath, project.downsample_targqc_dirpath)

        self.basecall_stat_html_reports = self.__get_basecall_stats_reports()
        info('basecall_stat_html_reports: ' + str(self.basecall_stat_html_reports))

        self.get_fastq_regexp_fn = get_hiseq4000_miseq_regexp

    def __find_fastq_dir(self):
        for dname in os.listdir(self.unaligned_dirpath):
            dpath = join(self.unaligned_dirpath, dname)
            if isdir(dpath) and any(f.endswith('.fastq.gz') for f in os.listdir(dpath)):
                return dpath

    def __get_basecall_stats_reports(self):
        self.basecalls_reports_dirpath = join(self.unaligned_dirpath, 'Reports', 'html')
        index_html_fpath = join(self.basecalls_reports_dirpath, 'index.html')
        if verify_dir(self.illumina_dir) and verify_file(index_html_fpath):
            return [index_html_fpath]


class NextSeq500Structure(DatasetStructure):
    def __init__(self, input_dir, proj_infos, samplesheet=None, **kwargs):
        info('Parsing the NextSeq500 project structure')
        self.kind = 'nextseq500'
        DatasetStructure.__init__(self, input_dir, proj_infos, samplesheet=samplesheet, **kwargs)
        debug('proj_infos: ' + str(proj_infos))

        verify_dir(self.unaligned_dirpath, is_critical=True)

        for pname, project in self.project_by_name.items():
            analysis_dir, output_dir, az_proj_name = self._get_output_dir(proj_infos, pname)

            project.set_dirpath(self.unaligned_dirpath, analysis_dir, output_dir, az_proj_name)
            for sample in project.sample_by_name.values():
                sample.source_fastq_dirpath = project.ds_dir
                sample.set_up_out_dirs(project.fastq_dirpath, project.fastqc_dirpath, project.downsample_targqc_dirpath)

        self.basecall_stat_html_reports = self.__get_basecall_stats_reports()

        self.get_fastq_regexp_fn = get_nextseq500_regexp

    def __find_fastq_dir(self):
        for dname in os.listdir(self.unaligned_dirpath):
            dpath = join(self.unaligned_dirpath, dname)
            if isdir(dpath) and any(f.endswith('.fastq.gz') for f in os.listdir(dpath)):
                return dpath

    def __get_basecall_stats_reports(self):
        self.basecalls_reports_dirpath = join(self.unaligned_dirpath, 'Reports', 'html')
        index_html_fpath = join(self.basecalls_reports_dirpath, 'index.html')
        if verify_dir(self.basecalls_reports_dirpath) and verify_file(index_html_fpath):
            return [index_html_fpath]


class DatasetProject:
    def __init__(self, name):
        self.name = name
        self.sample_by_name = OrderedDict()
        self.ds_dir = None
        self.output_dir = None
        self.analysis_dir = None
        self.az_project_name = None

        self.fastq_dirpath = None
        self.fastqc_dirpath = None
        self.downsample_metamapping_dirpath = None
        self.downsample_targqc_dirpath = None
        self.downsample_targqc_report_fpath = None
        self.multiqc_report_html_fpath = None
        self.mergred_dir_found = False

    def set_dirpath(self, ds_dir, analysis_dir, output_dir, az_project_name):
        self.ds_dir = ds_dir
        self.output_dir = output_dir.format(ds_proj_name=self.name)
        self.analysis_dir = analysis_dir
        self.az_project_name = az_project_name
        verify_dir(self.ds_dir, is_critical=True)

        found_merged_dirpath = join(self.ds_dir, 'merged')
        if verify_dir(found_merged_dirpath, silent=True):
            self.mergred_dir_found = True
            self.fastq_dirpath = self.fastqc_dirpath = found_merged_dirpath
        else:
            self.fastq_dirpath = join(self.ds_dir, 'fastq')
            try:
                safe_mkdir(self.fastq_dirpath)
            except:
                self.fastq_dirpath = join(self.output_dir, 'fastq')
            self.fastqc_dirpath = join(self.output_dir, 'FastQC')
        info()

        self.downsample_targqc_report_fpath = None
        self.multiqc_report_html_fpath = join(self.output_dir, 'multiqc_report.html')

        self.downsample_metamapping_dirpath = join(self.output_dir, 'Downsample_MetaMapping')
        self.downsample_targqc_dirpath = join(self.output_dir, 'Downsample_TargQC')
        self.downsample_targqc_report_fpath = join(self.downsample_targqc_dirpath, 'summary.html')

    def concat_fastqs(self, get_fastq_regexp):
        info('Preparing fastq files for the project named ' + self.name or self.az_project_name)
        if self.mergred_dir_found:
            info('  found already merged fastq dir, skipping.')
            return
        if not self.sample_by_name:
            err('  no samples found.')
            return
        try:
            if not isdir(self.fastq_dirpath):
                safe_mkdir(self.fastq_dirpath)
                open(join(self.fastq_dirpath, 'folder_created_by_prealign'), 'a').close()
        except OSError as e:
            err('Error: cannot write to ' + self.fastq_dirpath + '(' + str(e) + '). ' +
                'Writing fastq to ' + join(self.output_dir, 'fastq') + ' instead')
            self.fastq_dirpath = join(self.output_dir, 'fastq')

        try:
            fastqc_symlink = join(self.fastq_dirpath, 'FastQC')
            if not exists(fastqc_symlink):
                os.symlink(self.fastqc_dirpath, fastqc_symlink)
        except OSError:
            pass

        for s in self.sample_by_name.values():
            _concat_fastq(s.find_raw_fastq(get_fastq_regexp, 'R1'), s.l_fpath)
            _concat_fastq(s.find_raw_fastq(get_fastq_regexp, 'R2'), s.r_fpath)
        info()

class DatasetSample:
    def __init__(self, name, index=None, source_fastq_dirpath=None):
        self.name = re.sub(r'[\W_]', r'_', name)
        self.index = index

        self.fastq_dirpath = None
        self.source_fastq_dirpath = source_fastq_dirpath
        self.l_fpath = None
        self.r_fpath = None

        self.fastqc_dirpath = None
        self.l_fastqc_base_name = None
        self.r_fastqc_base_name = None
        self.l_fqc_sample = None
        self.r_fqc_sample = None

        self.lane_numbers = set()
        self.fcid = None  # for HiSeq

        self.targqc_sample = None
        self.downsample_targqc_dirpath = None

    def set_up_out_dirs(self, fastq_dirpath, fastqc_dirpath, downsample_targqc_dirpath):
        self.fastq_dirpath = fastq_dirpath
        self.fastqc_dirpath = fastqc_dirpath
        self.downsample_targqc_dirpath = downsample_targqc_dirpath

        self.l_fpath = join(fastq_dirpath, self.name + '_R1.fastq.gz')
        self.r_fpath = join(fastq_dirpath, self.name + '_R2.fastq.gz')

        # self.sample_fastqc_dirpath = join(fastqc_dirpath, self.name + '.fq_fastqc')
        # self.fastqc_html_fpath = join(fastqc_dirpath, self.name + '.fq_fastqc.html')
        self.l_fastqc_base_name = splitext_plus(basename(self.l_fpath))[0]
        self.r_fastqc_base_name = splitext_plus(basename(self.r_fpath))[0]
        # self.l_fastqc_html_fpath = None  # join(ds.fastqc_dirpath,  + '_fastqc.html')
        # self.r_fastqc_html_fpath = None  # join(ds.fastqc_dirpath, splitext_plus(self.r_fpath)[0] + '_fastqc.html')

        # if not isfile(self.fastqc_html_fpath):
        #     self.fastqc_html_fpath = join(self.sample_fastqc_dirpath, 'fastqc_report.html')

        # self.targqc_sample = targqc.Sample(self.name, join(downsample_targqc_dirpath, self.name), )
        # self.targqc_html_fpath = self.targqc_sample.targqc_html_fpath

    def find_raw_fastq(self, get_regexp, suf='R1'):
        fastq_fpaths = [
            join(self.source_fastq_dirpath, fname)
                for fname in os.listdir(self.source_fastq_dirpath)
                if re.match(get_regexp(self, suf), fname)]
        fastq_fpaths = sorted(fastq_fpaths)
        if not fastq_fpaths:
            critical('Error: no fastq files for the sample ' + self.name +
                     ' were found inside ' + self.source_fastq_dirpath)
        info(self.name + ': found raw fastq files ' + ', '.join(fastq_fpaths))
        return fastq_fpaths


def _concat_fastq(fastq_fpaths, output_fpath):
    if len(fastq_fpaths) == 1:
        if not isfile(output_fpath):
            info('  no need to merge - symlinking ' + fastq_fpaths[0] + ' -> ' + output_fpath)
            if not isdir(dirname(output_fpath)):
                safe_mkdir(output_fpath)
                critical('Dir for the symlink ' + dirname(output_fpath) + ' does not exist')
            os.symlink(fastq_fpaths[0], output_fpath)
            return output_fpath
    else:
        info('  merging ' + ', '.join(fastq_fpaths))
        if can_reuse(output_fpath, fastq_fpaths):
            info(output_fpath + ' exists, reusing')
        else:
            with file_transaction(None, output_fpath) as tx:
                with open(tx, 'wb') as out:
                    for fq_fpath in fastq_fpaths:
                        with open(fq_fpath, 'rb') as inp:
                            shutil.copyfileobj(inp, out)
        return output_fpath
