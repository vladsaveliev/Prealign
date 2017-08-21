## Pre-alignment processing and QC
This pipeline is ment to be run after the basecalling step. Make sure the bcl2fastq software was ran for the project (see BaseCalling / BaseCalling for HiSeq4000).

The following script detects the sequencing projects and samples from the SampleSheet.tsv file and the project structure (automatically detecting the type by base directory name: HiSeq / MiSeq / HiSeq4000). Then it does the following steps:

- Get the project name from the JIRA case, and create a symlink in the datasets folder as datasets/<division>/<project_name>
- If the fastq files are not yet merged by lanes, it merges them and saves into `<output_dir>/fastq`, or if output_dir is not specified: 
  - `Unalign/<project>/fastq` (for HiSeq4000)
  - `Unalign/Project_<project>/fastq` (for HiSeq)
  - `Data/Intensities/BaseCalls/fastq` (for MiSeq)
- Runs FastQC for each merged fastq file (separately for left and right reads), then combines them into one per project and saves into fastq/FastQC, and exposes both combined and separate htmls to the NGS webserver.
- Downsamples FastQC files to random 500.000 read pairs, aligns them to the reference, then runs and exposes the TargQC reports on resulting BAM files.
- Exposes the bcl2fastq reports:
  - `Unalign/Reports/html/index.html` (for MiSeq and HiSeq4000)
  - `Unalign/Basecall_Stats_<id>/All.htm`, `Unalign/Basecall_Stats_<id>/IVC.htm`, `Unalign/Basecall_Stats_<id>/Demultiplex_Stats.htm`, (for HiSeq).

##### Usage
```
module load preproc
prealign \
   /ngs/oncology/datasets/HiSeq/150612_D00443_0168_AHMNFGADXX \
   -o /ngs/oncology/analysis/dev/Dev_0104_HiSeq_DS \
   --genome hg19 \
   [--bed target.bed] \
   [--jira https://jira.rd.astrazeneca.net/browse/NGSG-313]
```

Arguments:
- First argument - path to the dataset (which contains 'Unalign' folder).
- `-o`: Desired output directory for merged reads and all reports, usually in "analysis".
- `--genome`: Genome version to align downsampled reads to (hg19, hg38, hg38-noalt)
- `--bed`: optional, used for downsampled TargQC reporting.
- `--jira`: JIRA case URL.

##### Multiple-project run
```
prealign \
   /ngs/oncology/Datasets/HiSeq4000/151222_K00172_0042_BH5GTFBBXX \
   --conf conf.csv
```

conf.csv
```
#project_name,output_dir,jira_url,bed
AZ300.Crown_Rerun,/ngs/oncology/analysis/dev/Dev_0138_HiSeq4000_CrownModels_Rerun,https://jira.rd.astrazeneca.net/i#browse/NGSG-440,target.bed
HyperPlus.Alpha_WGS,/ngs/oncology/analysis/dev/Dev_0165_HiSeq4000_AKT.WGS,https://jira.rd.astrazeneca.net/i#browse/NGSG-436,
```

Optionally, you can process a subproject separately by providing the exact subproject path in Unalign:
```
prealign \
   /ngs/oncology/datasets/HiSeq/150612_D00443_0168_AHMNFGADXX/Unalign/HyperPlus.Alpha_WGS \
   [--jira https://jira.rd.astrazeneca.net/i#browse/NGSG-436] \
   [--project-name Dev_0165_HiSeq4000_AKT.WGS3]
```

##### UK specific notes
Run bcl2fastq as 'sbsuser' first: http://wiki.rd.astrazeneca.net/display/NG/NextSeq+500+-+Basecalling+UK

As 'sbsuser', run
```
screen
module load preproc
preproc /ngs/oncology/datasets/NextSeq500/160318_NB501188_0011_AH3TCVAFXX --jira https://jira.rd.astrazeneca.net/browse/NGSG-536 --genome hg38 --project-name TS_UK_0010_Bladder_2_preqc
```
where the date stamped folder, Jira URL and project-name must be changed accordingly.
