#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="pymgal",
    version="1.0.0",
    author="Weiguang Cui",
    author_email="cuiweiguang@gmail.com",
    description="A Package for Mock Observations",
    long_description=read('README.md'),
    packages=find_packages(),
    requires=['numpy', 'pyfits', 'scipy', 'astropy'],
    package_data={
        '': ['*.fits',
             '*README*',
             'models/*.model',
             'filters/*',
             'refs/*']},
    license="BSD",
    include_package_data=True,
    keywords='astronomy astrophysics hydro-dynamical simulation mock observation',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3"
    ]
)
