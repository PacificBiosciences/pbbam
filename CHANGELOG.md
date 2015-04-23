# PacBio::BAM - change log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/). 

**NOTE:** The current series (0.y.z) is under initial development. Anything may change at any time. 
The public API should not be considered stable yet. Once we lock down a version 1.0.0, this will 
define a reference point & compatibility guarantees will be maintained within each major version 
series.


## [0.0.4] - 2015-04-22

### Added

- This changelog. Hope it helps.
- Hook to set verbosity of underlying htslib warnings.
- bug 26361 - grouped querys

### Changed

- Now using exceptions instead of return codes, output parameters, etc.
- Removed "messy" shared_ptrs across interface (see especially BamHeader). These are now taken care of within the API, not exposed to client code.

### Removed

- BamReader 

### Fixed

- bug 26381 - ASCII tag output
