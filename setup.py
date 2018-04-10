#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="sigseekr",
    version="0.2.0",
    packages=find_packages(),
    scripts=['sigseekr/sigseekr.py'],
    author="Andrew Low",
    author_email="andrew.low@inspection.gc.ca",
    url="https://github.com/lowandrew/SigSeekr",
    install_requires=['biopython', 'OLCTools', 'pytest']
)
