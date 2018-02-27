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

## FAQ

### [Help! I am getting "unsupported sequencing chemistry combination"!](#chemistry-bundle)

**pbbam** validates all BAM files, and as part of this validation, it checks whether the
`BindingKit` and `SequencingKit` variables in every ReadGroup of the provided BAM file are
known. As part of ongoing chemistry developments, we might need to introduce new part numbers
to identify novel reagents and/or SMRT Cells. You are unlikely to encounter such issues
when using SMRT Link, as it has an integrated auto-updater that will periodically check and
install new chemistries automatically. All PacBio tools being used without a proper SMRT Link
installation will require manual intervention to download new chemistries:

  ```sh
  cd <some persistent dir>
  export SMRT_CHEMISTRY_BUNDLE_DIR="${PWD}"

  wget https://raw.githubusercontent.com/PacificBiosciences/pbcore/develop/pbcore/chemistry/resources/mapping.xml -O chemistry.xml
  ```

This will cause **pbbam** to try to load the out-of-band `chemistry.xml` from
`SMRT_CHEMISTRY_BUNDLE_DIR` and should allow you to use somewhat older software
with somewhat newer BAMs. **Note:** this only allows **pbbam**'s internal validation
to pass, this will not automatically make other chemistry-dependent software work
with newer chemistries. For instance, Arrow's backend ([Unanimity](https://github.com/PacificBiosciences/unanimity))
is parametrized on chemistry too, and it will fail should a completely new chemistry
be introduced. See Unanimity's FAQ on how to employ `SMRT_CHEMISTRY_BUNDLE_DIR`
to load models for new chemistries.


## License

 - [PacBio open source license](https://github.com/PacificBiosciences/pbbam/blob/master/LICENSE.txt)

DISCLAIMER
----------
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

