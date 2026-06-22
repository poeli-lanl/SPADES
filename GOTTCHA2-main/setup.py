#!/usr/bin/env python3

from pathlib import Path
import re
from setuptools import find_packages, setup

ROOT = Path(__file__).parent.resolve()

def read_text(path):
    return (ROOT / path).read_text(encoding="utf-8")

def read_version():
    content = read_text("gottcha/gottcha2.py")
    match = re.search(r'^__version__\s*=\s*"([^"]+)"', content, re.MULTILINE)
    if not match:
        raise RuntimeError("Unable to find __version__ in gottcha/gottcha2.py")
    return match.group(1)

setup(
    name="gottcha2",
    version=read_version(),
    description="GOTTCHA2: Genomic Origin Through Taxonomic CHAllenge v2",
    long_description=read_text("README.md"),
    long_description_content_type="text/markdown",
    python_requires=">=3.9",
    author="Po-E Li",
    author_email="po-e@lanl.gov",
    license="BSD-3-Clause",
    url="https://github.com/poeli/GOTTCHA2",
    project_urls={
        "Homepage": "https://github.com/poeli/GOTTCHA2",
        "Source": "https://github.com/poeli/GOTTCHA2",
        "Issues": "https://github.com/poeli/GOTTCHA2/issues",
        "Documentation": "https://github.com/poeli/GOTTCHA2/blob/master/README.md",
    },
    keywords=[
        "bioinformatics",
        "taxonomy",
        "profiler",
        "metagenomics",
        "microbiome",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bioinformatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.2.0",
        "requests>=2.25.0",
        "biom-format>=2.1.7",
        "pysam>=0.22.0",
        "tqdm",
    ],
    scripts=['gottcha/utils/gottcha2.py', 'gottcha/utils/ont_utils.py', 'gottcha/utils/taxonomy.py'],
    packages=find_packages(include=["gottcha*"], exclude=["test*", "tests*", "docs*"]),
    include_package_data=True,
    zip_safe=False,
    entry_points={"console_scripts": ["gottcha2=gottcha.gottcha2:cli"]},
)
