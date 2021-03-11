# PacBio::BAM - change log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Active

### Fixed
 - Query name parsing for CCS/transcript records in PbiQueryNameFilter

## [1.6.2] - 2021-02-16

### Added
 - Use multi-threaded BAM reading,
   configurable via env variable PB_BAMREADER_THREADS

### Changed
 - Do not enforce subreads to have kinetics when converting BamRecord to Read

## [1.6.1] - 2020-11-10

### Added
 - Support for PBI index format version 4.0.0.
 - "Estimated bytes used" diagnostic for BamRecords. Memory usage is heavily
   implementation-dependent, w.r.t. data structure layout and alignment. This
   diagnostic is provided as a lower bound, but **no guarantee** is made for the
   exact value.

### Removed
 - Unused experimental 'CCSRecord' API

## [1.6.0] - 2020-09-23

### Added
 - ccs-kinetics-bystrandify tool.
 - 'errno' context for more informative I/O failure messages.
 - Support for alternative frame encodings.

### Changed
 - Moved BedReader/Writer to PacBio::BED namespace.

## [1.5.0] - 2020-08-05

### Added
 - FastaCache::Check methods for validating input FASTA
 - Detection of empty SAM/BAM input, as distinct from "could not read header"
   errors).
 - Full htslib v1.10 compatibility

### Fixed
 - Updating BAM record name no longer discards CCS strand suffixes.

### Changed
 - IndexedFastaReader no longer implicitly creates FASTA index (*.fai).

### Removed
 - 'HtslibVerbosity' setting in pbbam config. htslib's own logging to stderr is
   simply disabled at startup by default. In the rare case that client code needs
   to override this behavior, use 'hts_set_log_level()` directly.
 - Unused BamFile::FirstAlignmentOffset() method

## [1.4.0] - 2020-05-22

### Added
 - SamReader

### Changed
 - Data::SNR (from pbcopper) values are now float, not double.

### Removed
 - BamRecordBuilder

## [1.3.0] - 2020-04-24

### Added
 - FaiIndex::Create(), instead of needing samtools externally.

### Changed
 - Deprecated FrameEncodingType (not yet removed). Use FrameCodec instead.

## [1.2.0] - 2020-04-10

### Added
 - Support for older Linux kernels (<3.17) with our Dataset UUIDs.

### Fixed
 - BgzipWriter: default compression level.

## [1.1.0] - 2020-03-12

### Added
 - Support for barcode-labeled read groups.
 - "Run metadata" handling for dataset XML.

## [1.0.7] - 2019-10-10

### Added
 - CCSRecord API to work with the minimally required data for CCS
 - TextFileReader & TextFileWriter for generic line files (plain text or gzipped)
 - BedReader & BedWriter for BED format support
 - Override switch to allow de-novo generated DataSets to print verbatim (possibly
   relative) file paths, not only absolute paths.

### Removed
 - Unused built-in "convenience" queries: BarcodeQuery, QNameQuery, ReadAccuracy,
   SubreadLengthQuery. The same functionality can recreated via providing a
   PbiFilterQuery with a PbiFilter describing the corresponding criteria.

## [1.0.6] - 2019-06-14

### Added
 - IFastaWriter & IFastqWriter abstract base classes

## [1.0.5] - 2019-06-11

### Changed
 - BAM tag lookup improvements under the hood.

## [1.0.4] - 2019-06-07

### Added
 - General-purpose BgzipWriter
 - BgzipFastaWriter and BgzipFastqWriter
 - Read-only view to read indices passing a PbiFilter
 - IPD field to SimpleRead

## [1.0.3] - 2019-05-20

### Added
 - IndexedFastqReader for random access to FASTQ subregions

### Fixed
 - MappedSimpleRead clipping: on disjoint aligned/requested regions and on requests larger
   than available sequence.

## [1.0.2] - 2019-05-10

### Added
 - Range-for iteration on FastaReader & FastqReader

## [1.0.1] - 2019-05-09

### Added
 - SimpleRead & MappedSimpleRead for htslib-free processing.

### Fixed
 - Incorrect type displayed in SAM output (pure-text) for floating-point values.

## [1.0.0] - 2019-04-22

### Changed
 - C++14 is now a *hard* minimum.

### Removed
 - Headers emulating C++14 features for C++11.

### Fixed
 - Inconsistent whitelist/blacklist filters in DataSet XML.

## [0.25.0] - 2019-04-11

### Changed
 - Requires C++14 at minimum.

### Fixed
 - Reading BioSample(s) elements from DataSet XML.

## [0.24.0] - 2019-04-05

### Added
 - Built-in support for dataset elements: BioSample(s) & DNABarcode(s).
 - BaiIndexCache for reusing data from *.bai files(s).
 - Support in GenomicIntervalQuery for new BaiIndexCache.

## [0.23.1] - 2019-03-21

### Added
 - Streamable BamReader (via stdin).
 - Enabled range-for on BamReader, compatible with the other *Query inputs.

## [0.23.0] - 2019-03-11

### Added
 - PbiIndexCache and FastaCache for reusing file data
 - BaiIndexedBamReader and GenomicIntervalQuery can be constructed without
   initial interval.

## [0.22.0] - 2019-02-11

### Fixed
 - Handles zero-length reads for stitching ZMW reads.
 - Clipping to query on reverse-strand aligned reads.
 - Removed UB in dataset API.

### Added
 - "exciseFlankingInserts" option for clipping reads w.r.t reference.

## [0.21.0] - 2018-12-21

### Added
 - New local context flags: ADAPTER_BEFORE_BAD and ADAPTER_AFTER_BAD.

### Changed
 - Current PacBioBAM spec now 3.0.7.

### Removed
 - CMake has been removed completely.

## [0.20.0] - 2018-10-03

### Added
 - Support for (optionally) barcode-labeled read group IDs.

## [0.19.0] - 2018-09-11

### Added
 - TranscriptAlignmentSet to XML support

## [0.17.0] - 2018-03-18

### Added
- CompressionLevel/NumThreads parameter implementation to PbiBuilder.
- Dataset ctor to PbiFileQuery.
- TranscriptSet to XML support.
- Auto-enabled "permissive CIGAR mode" for pbbamify tool.
- IndexedBamWriter, for more efficient writing of BAM & PBI simultaneously.

## [0.16.0] - 2018-01-17

### Removed
- Removed the PbiIndex class and its "lookup data"-related helpers. These were
never as useful as initially intended. PbiRawData and its related classes are the
recommended interface for working with PBI index data.

## [0.15.0] - 2018-01-12

### Added
- Support for long CIGARs (>64K operations).

## [0.14.0] - 2017-12-12

### Added
- Support for newer style QNAMEs. Recent version of htslib (1.4+) have started
adding extra null terminators to make the subsequent CIGAR section 32-bit aligned.

### Changed
- Requirements for htslib version used. Must now be htslib v1.4+.

## [0.13.2] - 2017-09-25

### Added
- Backward compatibility for C++11 (std::make_unique which is 11/14 agnostic).

## [0.13.1] - 2017-09-25

### Added
- Support for "pe" tag in stitched, virtual reads.

## [0.13.0] - 2017-09-25

### Changed
- Ran clang-tidy (modernize) over codebase to clean up legacy coding styles.

## [0.12.2] - 2017-09-22

### Added
- HasPulseExclusion() to BamRecord (& derived types).

## [0.12.1] - 2017-09-21

### Added
- Pulse exclusion base feature to read group.

## [0.12.0] - 2017-09-19

### Added
- NumReads() for PBI filter-based queries. This allows fetching of the number
of reads that pass the filter, without needing to iterate over the entire
file(s).

## [0.11.0] - 2017-09-15

### Added
- Support for internal tag: pulse exclusion reason ("pe"). New methods on
BamRecord, and new enum PulseExclusionReason.

### Changed
- Default PacBioBAM format version now 3.0.5

## [0.10.2] - 2017-09-14

### Changed
- Explicitly trim all whitespace from FASTA input.

## [0.10.1] - 2017-09-11

### Changed
- Frames, add mutex to avoid race condition in InitIpdDownsampling(void)

## [0.10.0] - 2017-09-08

### Changed
- PbiBuilder backend for generating PBI index files "on-the-fly" along with
writing BAM files. The previous implementation's memory usage scaled linearly
with the number of reads, sometimes reaching huge numbers (several gigs or more).
The new implementation's memory usage remains constant for any number of reads,
without any runtime hit on files/architectures tested.

### Removed
- PbiBuilder::Result(). Returned an intermediate snapshot of the index under
construction. This method isn't usable with the new PbiBuilder backend and was
really only useful for initial debugging/testing. It is no longer used in the
test framework and is unlikely to be used by client code either. Dropping this
method from the API, and thus bumping the version number.

## [0.9.0] - 2017-08-07

### Removed
- Bundled htslib. Now using 'stock' htslib (v1.3.1+).
- Built-in SWIG wrappers.

## [0.8.0] - 2017-07-24

### Added
- Default DataSet 'Version' attribute if none already present (currently 4.0.0)
- Added whitelist support for filtering ZMWs via DataSetXML.
- Added iterable query over FASTA files & ReferenceSet datasets.
- Added DataSet::AllFiles to access primary resources AND their child files (indices,
scraps, etc).

### Fixed
- Bug in the build system preventing clean rebuilds.

### Removed
- Dropped the bundled, PacBio-forked version of htslib. Now using stock htslib (v1.3.1+).

## [0.7.4] - 2016-11-18

### Changed
- Compatibility for merging BAM files no longer requires exact match of PacBioBAM
version number (header @HD:pb tag). As long as both files meet the minimum
supported version number, the merge is allowed.

## [0.7.3] - 2016-11-11

### Added
- Support for S/P2-C2 chemistry and forthcoming 4.0 basecaller

## [0.7.2] - 2016-11-10

### Removed
- SAM header version equality check for merging BAM files. PacBioBAM version
number carries more meaning for PacBio data and thus will be the basis of
ensuring compatible merging.

## [0.7.1] - 2016-11-09

### Added
- (Unindexed) FASTA reader & FastaSequence data structure.
- Missing unit tests for internal BAM tag access.
- Chemistry data for basecaller v3.3.
- Missing parsers for filtering barcode quality ("bq"), barcode forward ("bcf"),
and barcode reverse ("bcr") from DataSetXML.
- Integrated htslib into project.

### Fixed
- Reverse complement on padding base.

## [0.7.0] - 2016-09-26

### Added
- Clipping for CCS records

### Fixed
- Cached position data leaking across records while iterating.
- Rolled back default pulse behavior in internal BAM API, to be backward-
compatible with existing client code (for now at least). v0.6.0 introduced
returning basecalled positions ONLY by default, rather than return ALL
pulses.
- Fixed crash when attempting to read from empty BAM/PBI files using the
PbiFilter-enabled APIs.

## [0.6.0] - 2016-09-13

### Added
- BamWriter writes to a BAM file with the target name plus a ".tmp" suffix. On
successful completion (i.e. normal BamWriter destruction, not triggered by a
thrown exception) the file is renamed to the actual requested filename.
- PBI file creation follows the same temporary naming convention.
- Support for barcode pair (forward, reverse) in DataSetXML filter.
- Validation API & 'auto-validate' compile-time switch.
- Added support for a batched QNAME whitelist filter in DataSet XML. Uses (new)
Property name 'qname_file', with the value being the filepath containing the
whitelist.
- Exposed MD5 hashing to API.
- Ability to remove base features from a ReadGroupInfo object.
- Can construct an aggregate PbiRawData index object from a DataSet: essentially
concatenates all PBI data within the dataset.
- New SamWriter class to create SAM-formatted output of PacBio BAM data.
- Extended APIs for accessing "internal BAM" data, including PulseBehavior
switch for selecting between all pulses & basecalls only.

### Fixed
- Improper 'clip to reference' product for BamRecord in some cases.
- Improper behavior in tag accessors (e.g. BamRecord::IPD()) on reverse strand-
aligned reads (bug 31339).
- Improper basecaller version parsing in ReadGroupInfo.

### Changed
- RecordType::POLYMERASE renamed to RecordType::ZMW to reflect changes in
PacBio BAM spec v3.0.4
- Refactored the 'virtual' reader classes - to match the new nomenclature,
and to combine the virtual reader & composite readers behind a shared
interface. The old class names still exist, as typedefs to the new ones,
and the interfaces are completely source-compatible - so as not to break
existing code. However, the old classes should be considered deprecated and
the new ones preferred. Below is the mapping of old -> new:

   VirtualPolymeraseBamRecord        ->  VirtualZmwBamRecord
   VirtualPolymeraseReader           ->  ZmwReadStitcher
   VirtualPolymeraseCompositeReader  ->  ZmwReadStitcher
   ZmwWhitelistVirtualReader         ->  WhitelistedZmwReadStitcher


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
