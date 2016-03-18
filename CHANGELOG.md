# PacBio::BAM - change log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/). 

**NOTE:** The current series (0.y.z) is under initial development. Anything may
change at any time. The public API should not be considered stable yet. Once we
lock down a version 1.0.0, this will define a reference point & compatibility
guarantees will be maintained within each major version series.

## Active

### Added
- BamWriter writes to a BAM file with the target name plus a ".tmp" suffix. On
successful completion (i.e. normal BamWriter destruction, not triggered by a
thrown exception) the file is renamed to the actual requested filename.
- PBI file creation follows the same temporary naming convention.
- Support for barcode pair (forward, reverse) in DataSetXML filter.
- Validation API & 'auto-validate' compile-time switch. 

### Fixed
- Improper 'clip to reference' product for BamRecord in some cases.

## [0.5.0] - 2016-02-22

### Added
- Platform model tag added to read group as RG::PM
- New scrap zmw type sz
- pbmerge accepts DataSetXML as input - using top-level resource BAMs as input,
applying filters, and generating a merged BAM. Also added FOFN support, instead
of listing out BAMs as command line args.
- PbiLocalContextFilter to allow filtering on subread local context.
- PbiBuilder: multithreading & zlib compression-level tuning for PBI output

### Fixed
- Fixed mishandling of relative BAM filenames in the filename constructor for
DataSet (e.g. DataSet ds("../data.bam")).

## [0.4.5] - 2016-01-14

### Changed
- PbiFilterQuery (and any other PBI-backed query, e.g. ZmwQuery ) now throws if
PBI file(s) missing insted of returning empty result.
- GenomicIntervalQuery now throws if BAI file(s) missing instead of returning
empty result.
- BamFile will throw if file is truncated (e.g. missing the EOF block). Disable
by defining PBBAM_NO_CHECK_EOF .

## [0.4.4] - 2016-01-07

### Added
- bam2sam command line utility. The primary benefit is removing the dependency
on samtools during tests, but also provides users a functioning BAM -> SAM
converter in the absence of samtools.
- pbmerge command line utility. Allows merging N BAM files into one, optionally
creating the PBI file alongside.
- Added BamRecord::Pkmean2 & Pkmid2, 2D equivalent of Pkmean/Pkmid, for internal
BAMs.

### Removed 
- samtools dependency

## [0.4.3] - 2015-12-22

### Added
- Compile using ccache by default, if available. Can be manually disabled using
-DPacBioBAM_use_ccache=OFF with cmake.
- pbindexdump: command-line utility that converts PBI file data into human-
readable formats. (JSON by default).

### Changed
- CMake option PacBioBAM_build_pbindex is being deprecated. Use
PacBioBAM_build_tools instead.

## [0.4.2] - 2015-12-22

### Changed
- BamFile::PacBioIndexExists & StandardIndexExists no longer check timestamps.
Copying/moving files around can yield timestamps that are not helpful (no longer
guaranteed that the .pbi will be "newer" than the .bam, even though no content
changed). Added methods (e.g. bool BamFile::PacBioIndexIsNewer()) to do that
lookup if needed, but it is no longer done automatically.

## [0.4.1] - 2015-12-18

### Added
- BamRecord::HasNumPasses

### Changed
- VirtualPolymeraseBamRecord::VirtualRegionsTable(type) returns an empty vector
of regions if none are associated with the requested type, instead of throwing.

## [0.4.0] - 2015-12-15

### Changed
- Redesigned PbiFilter interface and backend. Previous implementation did not
scale well as intermediate results were far too unwieldy. This redesign provides
speedups of orders of magnitude in many cases.

## [0.3.2] - 2015-12-10

### Added 
- Support for ReadGroupInfo sequencing chemistry data.
InvalidSequencingChemistryException thrown if an unsupported combination is
encountered.
- VirtualPolymeraseCompositeReader - for re-stitching records, across multiple
resources (e.g. from DataSetXML). Reader respects DataSet filter criteria.

## [0.3.1] - 2015-10-30

### Added
- ZmwWhitelistVirtualReader: similar to VirtualPolymeraseReader but restricts
iteration to a whitelist of ZMW hole numbers, leveraging PBI index data for
random-access.

### Fixed
- Fixed error in PBI construction, in which entire file sections (e.g.
BarcodeData or MappedData) where being dropped when any one record lacked data.
Correct behavior is to allow file section ommission if all records lack that
data type.

## [0.3.0] - 2015-10-29

### Fixed
- Improper reporting of current offset from multi-threaded BamWriter. This had
the effect of creating broken PBIs that were written alongside the BAM. Added a
flush step, which incurs a performance hit, but restores correctness.

## [0.2.4] - 2015-10-26

### Fixed
- Empty PbiFilter now returns all records, instead of filtering away all records.

## [0.2.3] - 2015-10-26

### Added/Fixed
- Syncing DataSetXML across APIs. Primary changes include output of Version
attribute ("3.0.1") on appropriate elements, as well as resolution of namespace
issues.

## [0.2.2] - 2015-10-22

### Added
- Added BAI bin calculation to BamWriter::Write, to ensure maximal compatibility
with downstream tools (e.g. 'samtools index'). A new BinCalculationMode enum
flag in BamWriter constructor cotnrols whether this behavior is enabled[default]
or not.

## [0.2.1] - 2015-10-19

### Added
- Exposed the following classes to public API:
  - BamReader
  - BaiIndexedBamReader
  - PbiIndexedBamReader
  - GenomicIntervalCompositeBamReader
  - PbiFilterCompositeBamReader

## [0.2.0] - 2015-10-09

### Changed
- BAM spec v3.0.1 compliance. Previous (betas) versions of the BAM spec are not
supported and will causean exception to be throw if encountered.
- PBI lookup interface & backend, see PbiIndex.h & PbiLookupData.h for details.

### Added 
- BamFile::PacBioIndexExists() & BamFile::StandardIndexExists() - query the
existence of index files without auto-building them if they are missing, as in
BamFile::Ensure*IndexExists().
- GenomicInterval now accepts an htslib/samtools-style REGION string in the
constructor: GenomicInterval("chr1:1000-2000"). Please note though, that pbbam
uses 0-based coordinates throughout, whereas samtools expects 1-based. The above
string is equivalent to "chr1:1001-2000" in samtools.
- Built-in PBI filters. See PbiFlter.h & PbiFilterTypes.h for built-in filters
and constructing composite filters. These can be used in conjunction with the
new PbiFilterQuery, which takes a generic PbiFilter and applies that to a
DataSet for iteration.
- New built-in queries: BarcodeQuery, ReadAccuracyQuery, SubreadLengthQuery.
These leverage the new filter API to construct a PbiFilter and apply to a
DataSet.
- Built-in BamRecord comparators that are STL-compatible. See Compare.h for full
list. This allows for statements like the following, which sorts records by ZMW
number:
``` c++
    vector<BamRecord> data;
    std::sort(data.begin(), data.end(), Compare::Zmw());
```
- "exciseSoftClips" option to BamRecord::CigarData()

## [0.1.0] - 2015-07-17

### Changed
- BAM spec v3.0b7 compliance
 - Removal of 'M' as allowed CIGAR operation. Attempt to use such a CIGAR op
 will throw an exception.
 - Addition of IPD/PulseWidth codec version info in header
  
### Added
- Auto-generation of UTC timestamp for DataSet objects
- PbiBuilder - allows generation of PBI index data alongside generation or
modification of BAM record data. This obviates the need to wait for a completed
BAM, then go through the zlib decompression, etc.
- Added DataSet::FromXml(string xml) to create DataSets from "raw" XML string,
rather than building up using DataSet API or loading from existing file.
- "pbindex" command line tool to generate ".pbi" files from BAM data. The
executable is built by default, but can be disabled using the cmake option
"-DPacBioBAM_build_pbindex=OFF".
  
### Fixed
- PBI construction failing on CCS reads

## [0.0.8] - 2015-07-02

### Changed
- Build system refactoring.

## [0.0.7] - 2015-07-02

### Added
- PBI index lookup API. Not so much intended for client use directly, but will
enable construction of higher-level semantic queries: grouping by, filtering,
etc.
- DataSet & PBI-aware queries (e.g. ZmwGroupQuery). More PBI-enabled queries to
follow.
- More flexibility in tag access. Samtools has a habit of performing a
"shrink-to-fit" when it handles integer-valued tag data. Thus we cannot
**guarantee** the binary type that our API will have to process. Safe
conversions are allowed on integer-like data only. Under- or overflows in
casting will trigger an exception. All other tag data types must be asked for
explicitly, or else an exception will be raised, as before.
- BamHeader::DeepCopy - allows creation of editable header data, without
overwriting all shared instances

### Fixed
- XSD compliance for DataSet APIs.

### Changed
- The functionality provided by ZmwQuery (group by hole number), is now
available using the ZmwGroupQuery object. The new ZmwQuery returns a single-
record iterator (a la EntireFileQuery), but limited to a whitelist of requested
hole numbers.

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

- DataSet support. This includes XML I/O, basic dataset query/manipulation, and
multi-BAM-file queries. New classes are located in <pbbam/dataset/>. DataSet-
capable queries currently reside in the PacBio::BAM::staging namespace. These
will be ported over to the main namespace once the support is stabilized and
works seamlessly with either a single BamFile or DataSet object as input. (bug
25941)
- PBI support. This includes read/write raw data & building from a BamFile. The
lookup API for random-access queries is under development, but the raw data is
available - for creating PBI files & generating summary statistics. (bug 26025)
- C# SWIG bindings, alongside existing Python and R wrappers.
- LocalContextFlags support in BamRecord (bug 26623)

### Fixed

- BamRecord[Impl] map quality now  initialized with 255 (missing) value, instead
of 0. (bug 26228)
- ReadGroupId calculation. (bug 25940)
  
## [0.0.4] - 2015-04-22

### Added

- This changelog. Hope it helps.
- Hook to set verbosity of underlying htslib warnings.
- Grouped queries. (bug 26361)

### Changed

- Now using exceptions instead of return codes, output parameters, etc.
- Removed "messy" shared_ptrs across interface (see especially BamHeader). These
are now taken care of within the API, not exposed to client code.

### Removed

- BamReader 

### Fixed

- ASCII tag output. (bug 26381)
