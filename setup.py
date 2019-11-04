#!/usr/bin/env python
import setuptools

with open("README.org", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="jrbp",
    version="0.0.1",
    author="John Bonini",
    author_email="john@rbonini.science",
    description="Python interface to ISOTROPY for Linux",
    long_description=long_description,
    long_description_content_type="text/org",
    url="https://github.com/jrbp/pysotropy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.6',
)
