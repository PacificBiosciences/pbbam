# pbbam

[![Build Status](https://travis-ci.org/PacificBiosciences/pbbam.svg?branch=master)](https://travis-ci.org/PacificBiosciences/pbbam) [![Documentation Status](https://readthedocs.org/projects/pbbam/badge/?version=latest)](http://pbbam.readthedocs.org/en/latest/?badge=latest)
 
As of the 3.0 release of SMRTanalysis, PacBio is embracing the industry standard BAM
format for (both aligned and unaligned) basecall data files. We have also formulated
a BAM companion file format (bam.pbi) enabling fast access to a richer set of per-read
information as well as compatibility for software built around the legacy cmp.h5 format.
 
The **pbbam** software package provides components to create, query, & edit PacBio BAM
files and associated indices. These components include a core C++ library, bindings for
additional languages, and command-line utilities.

### Note:

This library is **not** intended to be used as a general-purpose BAM utility - all input & output BAMs must adhere to the [PacBio BAM format specification](https://github.com/PacificBiosciences/PacBioFileFormats/blob/3.0/BAM.rst). Non-PacBio BAMs will cause exceptions to be thrown.
 
##  Documentation

  - [Documentation Home](http://pbbam.readthedocs.org/en/latest/index.html)
    - [Getting Started](http://pbbam.readthedocs.org/en/latest/getting_started.html)
    - [C++ API Reference](http://pbbam.readthedocs.org/en/latest/api_reference.html)

  - [Changelog](https://github.com/PacificBiosciences/pbbam/blob/master/CHANGELOG.md)

## License

 - [PacBio open source license](https://github.com/PacificBiosciences/pbbam/blob/master/LICENSE.txt)

DISCLAIMER
----------
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

