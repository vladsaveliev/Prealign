#!/usr/bin/env python
import sys
py_v = sys.version_info[:2]
if not (py_v == (2, 7) or py_v >= (3, 3)):
    sys.exit('Only Python 2.7 or 3.3 and up are supported. Current version: ' + '.'.join(py_v))

import os
import sys

package_name = script_name = 'prealign'


from setuptools import setup

try:
    from ngs_utils import setup_utils
except:
    version = open('VERSION.txt').read().strip().split('\n')[0]
    setup(version=version)  # For conda-build jinja context to read the version
else:
    version = setup_utils.init(package_name, package_name, __file__)
    setup(
        name=package_name,
        version=version,
        author='Vlad Saveliev and Alla Mikheenko',
        author_email='vladislav.sav@gmail.com',
        description='AstraZeneca pre-alignment analysis and reporting suite',
        long_description=(open('README.md').read()),
        keywords='bioinformatics',
        url='https://github.com/AstraZeneca-NGS/Pre_Alignment_Suite',
        download_url='https://github.com/AstraZeneca-NGS/Pre_Alignment_Suite/releases',
        license='GPLv3',
        packages=[package_name],
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
        scripts=[os.path.join('scripts', script_name)],
        include_package_data=True,
        zip_safe=False,
        install_requires=setup_utils.get_reqs(),
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
