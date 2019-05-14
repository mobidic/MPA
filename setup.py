#!/usr/bin/env python3
#
# Copyright (C) 2019
#

__author__ = 'Mobidic'
__authors__ = ['Henri Pegeot','Kevin Yauy','Charles Van Goethem','David Baux']
__copyright__ = 'Copyright (C) 2019'
__license__ = 'Academic License Agreement'
__version__ = '1.0.0'
__email__ = 'c-vangoethem@chu-montpellier.fr'
__status__ = 'dev'

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mobidic-mpa",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="MPA: MoBiDiC Prioritization Algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://neuro-2.iurc.montp.inserm.fr/mpaweb/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=['pyvcf==0.6.8'],
    entry_points={
        "console_scripts": [
            "mpa_main=mobidic_mpa_test:main"
        ],
    },
    scripts = ['scripts/mpa'],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/mobidic/MPA/issues',
        'Source': 'https://github.com/mobidic/MPA',
},
)
