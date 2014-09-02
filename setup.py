from ez_setup import use_setuptools
use_setuptools()

import re

from setuptools import setup, find_packages


def version_handler(mgr, options):
    version = mgr.get_current_version()
    if version.endswith('dev'):
        version += '-' + mgr._invoke('log', '-l1', '-r.', '--template', '{node|short}').strip()
    elif re.match('^\d+\.\d+$', version):
        # StrictVersion considers x.y == x.y.0 and drops the .0 from a
        # repo tag.  Add it back and ensure that it's really a tag for
        # our parent.
        version += '.0'
        assert version in mgr.get_parent_tags('tip')
    return version


setup(
    license = 'Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)',
    long_description=open('README.rst').read(),
    use_vcs_version = {'version_handler': version_handler},
    zip_safe = True,

    author = 'Reece Hart',
    author_email = 'reecehart@gmail.com',
    description = """miscellaneous simple bioinformatics tools and lookup tables""",
    name = "biotools",
    packages = find_packages(),
    url = 'https://bitbucket.org/reece/bioutils',

    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python",
        "Topic :: Database :: Front-Ends",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],

    keywords = [
        'bioinformatics',
    ],

    install_requires = [
        'biopython',
        'pysam >= 0.8.0',
    ],

    setup_requires = [
        'hgtools',
        'nose',
        #'sphinx',
        #'sphinxcontrib-fulltoc',
    ],

    tests_require = [
        #'coverage',
    ],
)

## <LICENSE>
## Copyright 2014 Bioutils Contributors (https://bitbucket.org/uta/bioutils)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
