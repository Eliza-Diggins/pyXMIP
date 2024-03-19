#!/usr/bin/env python
from setuptools import setup

setup(
    name="pyXs",
    packages=["pyXs"],
    version="0.1.0",
    description="",
    author="Eliza C. Diggins",
    author_email="eliza.diggins@utah.edu",
    url="https://github.com/Wik-Group/pyXs",
    download_url="https://github.com/Wik-Group/pyXs/tarball/0.1.0",
    install_requires=["numpy", "scipy", "astropy", "astroquery", "rich"],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    include_package_data=True,
)
