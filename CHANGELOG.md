# PacBio::BAM - change log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/). 

**NOTE:** The current series (0.y.z) is under initial development. Anything may change at any time. 
The public API should not be considered stable yet. Once we lock down a version 1.0.0, this will 
define a reference point & compatibility guarantees will be maintained within each major version 
series.

## [0.0.5] - 2015-05-29

### Added

- DataSet support. This includes XML I/O, basic dataset query/manipulation, and multi-BAM-file 
  queries. New classes are located in <pbbam/dataset/>. DataSet-capable queries currently reside in the 
  PacBio::BAM::staging namespace. These will be ported over to the main namespace once the support is 
  stabilized and works seamlessly with either a single BamFile or DataSet object as input. (bug 25941)
- PBI support. This includes read/write raw data & building from a BamFile. The lookup API for 
  random-access queries is under development, but the raw data is available - for creating PBI files & 
  generating summary statistics. (bug 26025)
- C# SWIG bindings, alongside existing Python and R wrappers.
- LocalContextFlags support in BamRecord (bug 26623)

### Fixed

- BamRecord[Impl] map quality now  initialized with 255 (missing) value, instead of 0. (bug 26228)
- ReadGroupId calculation. (bug 25940)
  
## [0.0.4] - 2015-04-22

### Added

- This changelog. Hope it helps.
- Hook to set verbosity of underlying htslib warnings.
- Grouped queries. (bug 26361)

### Changed

- Now using exceptions instead of return codes, output parameters, etc.
- Removed "messy" shared_ptrs across interface (see especially BamHeader). These are now taken care of within the API, not exposed to client code.

### Removed

- BamReader 

### Fixed

- ASCII tag output. (bug 26381)
