#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
from setuptools import setup, find_packages, find_namespace_packages
try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements


MAJOR = 1
MINOR = 0
VERSION = '%d.%d' % (MAJOR, MINOR)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def write_version_py(filename='pymgal/__version__.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SETUP.PY
__version__ = '%(version)s'
git_revision = '%(git_revision)s'
"""
    GIT_REVISION = git_version()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'git_revision': GIT_REVISION})
    finally:
        a.close()


def setup_package():
    with open('requirements.txt') as f:
        requirements = f.read().splitlines()
    write_version_py()
    setup(
        name='pymgal',
        packages=find_namespace_packages(),
        version=VERSION,
        author='Weiguang Cui, Patrick Janulewicz',
        author_email='cuiweiguang@gmail.com, patrick.janulewicz@gmail.com',
        description="A Package for Mock Observations in optical bands",
        long_description=read('README.md'),
        long_description_content_type='text/markdown',
        #packages=find_namespace_packages(),
        # requires=['numpy', 'scipy', 'astropy', 'h5py', 'PyYAML'],
        install_requires = requirements, 
        #parse_requirements('requirements.txt', session='hack'),
        #package_dir={"./": "./pymgal/"},
        #package_data={"models": ["*.model"],"filters": ["*"]},
        include_package_data=True,
        package_data={
            '': ['*.fits',
                 '*README*',
                 'models/*.model',
                 'filters/*',
                 'refs/*']},
        license="BSD",
        keywords='astronomy astrophysics hydrodynamical simulation mock observation',
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


if __name__ == '__main__':
    setup_package()