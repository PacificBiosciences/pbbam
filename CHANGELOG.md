# PacBio::BAM - change log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/). 

**NOTE:** The current series (0.y.z) is under initial development. Anything may change at any time. 
The public API should not be considered stable yet. Once we lock down a version 1.0.0, this will 
define a reference point & compatibility guarantees will be maintained within each major version 
series.

## Active

## [0.2.0] - 2015-10-09

### Changed
- BAM spec v3.0.1 compliance. Previous (betas) versions of the BAM spec are not supported and will cause
  an exception to be throw if encountered.
- PBI lookup interface & backend, see PbiIndex.h & PbiLookupData.h for details.

### Added 
- BamFile::PacBioIndexExists() & BamFile::StandardIndexExists() - query the existence of index files 
without auto-building them if they are missing, as in BamFile::Ensure*IndexExists().
- GenomicInterval now accepts an htslib/samtools-style REGION string in the constructor: 
GenomicInterval("chr1:1000-2000"). Please note though, that pbbam uses 0-based coordinates throughout, 
whereas samtools expects 1-based. The above string is equivalent to "chr1:1001-2000" in samtools.
- Built-in PBI filters. See PbiFlter.h & PbiFilterTypes.h for built-in filters and constructing composite
  filters. These can be used in conjunction with the new PbiFilterQuery, which takes a generic PbiFilter
  and applies that to a DataSet for iteration.
- New built-in queries: BarcodeQuery, ReadAccuracyQuery, SubreadLengthQuery. These leverage the new filter
  API to construct a PbiFilter and apply to a DataSet.
- Built-in BamRecord comparators that are STL-compatible. See Compare.h for full list. 
  This allows statements like:
      vector<BamRecord> data;
      std::sort(data.begin(), data.end(), Compare::Zmw());
  to get a container of records, sorted by ZMW. 
- "exciseSoftClips" option to BamRecord::CigarData() 

## [0.1.0] - 2015-07-17

### Changed
- BAM spec v3.0b7 compliance
 - Removal of 'M' as allowed CIGAR operation. Attempt to use such a CIGAR op will throw an exception.
 - Addition of IPD/PulseWidth codec version info in header
  
### Added
- Auto-generation of UTC timestamp for DataSet objects
- PbiBuilder - allows generation of PBI index data alongside generation/modification of BAM record
data. This obviates the need to wait for a completed BAM, then go through the zlib decompression, etc.
- Added DataSet::FromXml(string xml) to create DataSets from "raw" XML string, rather than building up 
using DataSet API or loading from existing file.
- "pbindex" command line tool to generate ".pbi" files from BAM data. The executable is built by default, 
but can be disabled using the cmake option "-DPacBioBAM_build_pbindex=OFF".
  
### Fixed
- PBI construction failing on CCS reads

## [0.0.8] - 2015-07-02

### Changed
- Build system refactoring.

## [0.0.7] - 2015-07-02

### Added
- PBI index lookup API. Not so much intended for client use directly, but will enable construction of
  higher-level semantic queries: grouping by, filtering, etc.
- DataSet & PBI-aware queries (e.g. ZmwGroupQuery). More PBI-enabled queries to follow.
- More flexibility in tag access. Samtools has a habit of performing a "shrink-to-fit" when it handles
  integer-valued tag data. Thus we cannot **guarantee** the binary type that our API will have to process.
  Safe conversions are allowed on integer-like data only. Under- or overflows in casting will trigger an 
  exception. All other tag data types must be asked for explicitly, or else an exception will be raised, 
  as before.
- BamHeader::DeepCopy - allows creation of editable header data, without overwriting all shared instances

### Fixed
- XSD compliance for DataSet APIs.

### Changed
- The functionality provided by ZmwQuery (group by hole number), is now available using the ZmwGroupQuery
  object. The new ZmwQuery returns a single-record iterator (a la EntireFileQuery), but limited to a whitelist 
  of requested hole numbers.

### Removed
- XSD non-compliant classes (e.g. ExternalDataReference)

## [0.0.6] - 2015-06-07

### Added

- Accessor methods for pulse bam support:
 - LabelQV()
 - AltLabelQV()
 - LabelTag()
 - AltLabelTag()
 - Pkmean()
 - Pkmid()
 - PrePulseFrames() only RC, no clipping
 - PulseCallWidth() only RC, no clipping
 - PulseCall() case-sensitive RC, no clipping
 - IPDRaw() to avoid up and downscaling for stitching
- BamRecord::ParseTagName and BamRecord::ParseTagString to convert a two 
  character tag string to a TagName enum and back. Allows a switch over tags.
- VirtualPolymeraseReader to create VirtualPolymeraseBamRecord from a 
  subreads|hqregion+scraps.bam
- VirtualRegion represents annotations of the polymerase reads, for adapters, 
  barcodes, lqregions, and hqregions.
- ReadGroupInfo operator== 

### Fixed

- Reimplemented QueryStart(int), QueryEnd(int), UpdateName(void), 
  ReadGroup(ReadGroupInfo&), ReadGroupId(std::string&);

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
