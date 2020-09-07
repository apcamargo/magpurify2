# !/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the magpurify2 package, available at:
# https://github.com/apcamargo/magpurify2
#
# Magpurify2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact: antoniop.camargo@gmail.com

from setuptools import find_packages, setup
from setuptools_rust import RustExtension

setup(
    name="magpurify2",
    version="0.1.0",
    packages=find_packages(),
    rust_extensions=[
        RustExtension("magpurify2._codon", debug=False),
        RustExtension("magpurify2._composition", debug=False),
        RustExtension("magpurify2._coverage", debug=False),
        RustExtension("magpurify2._mmseqs2", debug=False),
    ],
    zip_safe=False,
    license="GNU General Public License v3.0",
    description="Identify and remove contaminants from metagenome-assembled genomes.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "biolib",
        "biopython",
        "hdbscan",
        "importlib-metadata>=0.12; python_version<'3.8'",
        "joblib",
        "numpy",
        "scipy",
        "taxopy",
        "tbb",
        "umap-learn",
        "xgboost",
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["magpurify2=magpurify2.cli:cli"]},
    url="https://github.com/apcamargo/magpurify2/",
    keywords=["bioinformatics", "metagenomics", "metagenome-assembled genomes"],
    author="Antonio Camargo",
    author_email="antoniop.camargo@gmail.com",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python :: 3",
    ],
)
