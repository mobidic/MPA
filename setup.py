#!/usr/bin/env python3.6
#
# Copyright (C) 2019
#

__author__ = ['Mobidic']
__authors__ = ['Henri Pegeot','Kevin Yauy','Charles Van Goethem']
__copyright__ = 'Copyright (C) 2019'
__license__ = 'Academic License Agreement'
__version__ = '0.3.0a0'
__email__ = 'h-pegeot@chu-montpellier.fr'
__status__ = 'prod'

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mobidic_mpa_1",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="MoBiDiC Prioritization Algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mobidic/MPA",
    packages=setuptools.find_packages(exclude=['contrib', 'docs', 'tests']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ["pyvcf"]
)
