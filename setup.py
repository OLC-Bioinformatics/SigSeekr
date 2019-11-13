#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="sigseekr",
    version="0.2.1",
    packages=find_packages(),
    scripts=['sigseekr/sigseekr.py'],
    author="Adam Koziol",
    author_email="adam.koziol@canada.ca",
    url="https://github.com/OLC-Bioinformatics/SigSeekr.git",
    install_requires=['biopython', 'OLCTools', 'pytest']
)
