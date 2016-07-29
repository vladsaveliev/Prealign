#!/usr/bin/env python

import os
import sys
from setuptools import setup, find_packages
from pip.req import parse_requirements

try:
    from pip.download import PipSession
except ImportError:  # newer setuptools
    install_reqs = parse_requirements("requirements.txt")
else:
    install_reqs = parse_requirements("requirements.txt", session=PipSession())

reqs = [str(ir.req) for ir in install_reqs]
version = open('VERSION.txt').read().strip()

name = 'AZ_PreAlignmentSuite'
script = 'prealign'

def _run(_cmd):
    print('$ ' + _cmd)
    os.system(_cmd)

if sys.argv[-1] == 'publish':
    _run('python setup.py sdist upload')
    _run('python setup.py bdist_wheel upload')
    sys.exit()

if sys.argv[-1] == 'tag':
    _run('git tag -a %s -m "Version %s"' % (version, version))
    _run("git push --tags")
    sys.exit()

if sys.argv[-1] == 'install':
    print("""-----------------------------------
 Installing AZ pre-alignment suite version {}
-----------------------------------
""".format(version))

setup(
    name=name,
    version=version,
    author='Vlad Saveliev and Alla Mikheenko',
    author_email='vladislav.sav@gmail.com',
    description='AstraZeneca pre-alignment analysis and reporting suite',
    long_description=(open('README.md').read()),
    keywords='bioinformatics',
    url='https://github.com/AstraZeneca-NGS/Pre_Alignment_Suite',
    download_url='https://github.com/AstraZeneca-NGS/Pre_Alignment_Suite/releases',
    license='GPLv3',
    packages=find_packages(),
    package_data={
        'Utils': [
            'reference_data/fai/*.fai',
            'reporting/static/*.js',
            'reporting/static/*/*.js',
            'reporting/static/*.css',
            'reporting/static/*/*.css',
            'reporting/static/*.json',
            'reporting/static/*/*.json',
            'reporting/static/*.png',
            'reporting/static/*/*.png',
            'reporting/static/*.pxm',
            'reporting/static/*/*.pxm',
            'reporting/*.html',
            'reporting/*.json',
            'sambamba_binaries/sambamba_*',
            'sambamba/build/sambamba',
            'tools/*.sh',
        ],
        'prealign': [
            'configs/*.yaml',
            'webserver/id_rsa',
            'webserver/id_rsa.pub',
        ]
    },
    include_package_data=True,
    zip_safe=False,
    scripts=['scripts/' + script],
    install_requires=reqs,
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)

if sys.argv[-1] == 'install':
    print("""
--------------------------------
 AZ pre-alignment suite installation complete!
--------------------------------
Usage: {script} \
   /ngs/oncology/datasets/HiSeq/150612_D00443_0168_AHMNFGADXX \
   [--jira https://jira.rd.astrazeneca.net/browse/NGSG-313] \
   [--bed target.bed] \
   [--project-name Dev_0104_HiSeq_DS]

For help in running the suite, please see the documentation available at https://github.com/AstraZeneca-NGS/Pre_Alignment_Suite or run: {script} --help
""".format(script=script))
