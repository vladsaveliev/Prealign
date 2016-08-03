import socket
from contextlib import contextmanager
from os.path import abspath, expanduser, join, dirname, pardir
from traceback import format_exc

from Utils.config import load_yaml_config
from Utils.logger import info, err, critical, debug, warn
from Utils.file_utils import verify_file, verify_module, adjust_path
from Utils.utils import is_uk, is_us, is_cloud, is_sweden, is_ace, is_china, is_local

from yaml import load as load_yaml
try:
    from yaml import CDumper as Dumper, CLoader as Loader
except ImportError:
    from yaml import Dumper, Loader

code_base_path = abspath(join(dirname(abspath(__file__))))
configs_dirpath = abspath(join(code_base_path, 'configs'))
test_dirpath = abspath(join(code_base_path, 'test'))
sys_cnfs = dict(
    us = join(configs_dirpath, 'system_info_Waltham.yaml'),
    uk = join(configs_dirpath, 'system_info_AP.yaml'),
    sweden = join(configs_dirpath, 'system_info_Sweden.yaml'),
    china = join(configs_dirpath, 'system_info_China.yaml'),
    cloud = join(configs_dirpath, 'system_info_cloud.yaml'),
    ace = join(configs_dirpath, 'system_info_ace.yaml'),
    local = join(configs_dirpath, 'system_info_local.yaml'),
)

downsample_pairs_num = 5e5
genome = 'hg19'
fai_fpath = None
dedup = True

reuse_intermediate = False
is_debug = False
parallel_cfg = None
sys_cfg = None


project_kind_by_prefix = {
    'Bio_': 'bioscience',
    'Dev_': 'dev',
    'EXT_': 'external',
    'TS_': 'translation',
}


def load_sys_cfg(genome_build):
    sys_cfg_fpath = verify_file(_detect_sys_cnf_by_location(), is_critical=True)
    sys_cfg_ = load_yaml_config(sys_cfg_fpath)
    if genome_build not in sys_cfg_['genomes']:
        critical(genome_build + ' is not supported')
    sys_cfg_['genome'] = sys_cfg_['genomes'][genome_build]
    sys_cfg_['genome']['name'] = genome_build
    return sys_cfg_


def _detect_sys_cnf_by_location():
    if is_uk():
        res = sys_cnfs['uk']
    elif is_sweden():
        res = sys_cnfs['sweden']
    elif is_china():
        res = sys_cnfs['china']
    elif is_us():
        res = sys_cnfs['us']
    elif is_cloud():
        res = sys_cnfs['cloud']
    elif is_local():
        res = sys_cnfs['local']
    elif is_ace():
        res = sys_cnfs['ace']
    else:
        warn('Warning: could not detect location by hostname: ' + socket.gethostname() + '. Using local')
        res = sys_cnfs['local']
    return res


@contextmanager
def with_cnf(cnf, **kwargs):
    prev_opts = {k: cnf[k] for k in kwargs.keys()}
    try:
        for k, v in kwargs.items():
            cnf[k] = v
        yield cnf
    finally:
        for k, v in prev_opts.items():
            cnf[k] = v