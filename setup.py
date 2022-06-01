#!/usr/bin/env python3
#
# Copyright (C) 2017-2022
#

__author__ = 'Mobidic'
__authors__ = [
    'Henri Pegeot',
    'Kevin Yauy',
    'Charles Van Goethem',
    'Thomas Guignard',
    'David Baux'
]
__copyright__ = 'Copyright (C) 2017-2022'
__license__ = 'Academic License Agreement'
__version__ = '1.3.0'
__email__ = 'c-vangoethem@chu-montpellier.fr'
__status__ = 'prod'

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
    url="https://mobidic.github.io/MPA/",
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
    install_requires=[
        'vcfpy>=0.13.4',
        'tqdm>=4.59.0',
        'pysam>=0.19.1'
    ],
    scripts=['scripts/mpa'],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/mobidic/MPA/issues',
        'Source': 'https://github.com/mobidic/MPA',
    }
)
