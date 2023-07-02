#ifndef PBBAM_BAMRECORD_H
#define PBBAM_BAMRECORD_H

#include <pbbam/Config.h>

#include <pbbam/BamHeader.h>
#include <pbbam/BamRecordImpl.h>
#include <pbbam/ClipType.h>
#include <pbbam/FrameEncodingType.h>
#include <pbbam/PulseBehavior.h>
#include <pbbam/PulseExclusionReason.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/RecordType.h>
#include <pbbam/ZmwType.h>
#include <pbbam/virtual/VirtualRegionType.h>

#include <pbcopper/data/Accuracy.h>
#include <pbcopper/data/Frames.h>
#include <pbcopper/data/LocalContextFlags.h>
#include <pbcopper/data/MappedRead.h>
#include <pbcopper/data/Orientation.h>
#include <pbcopper/data/QualityValues.h>
#include <pbcopper/data/Read.h>
#include <pbcopper/data/Strand.h>
#include <pbcopper/json/JSON.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

class Pulse2BaseCache;

/// \brief The BamRecord class represents a %PacBio %BAM record.
///
/// %PacBio %BAM records are extensions of normal SAM/BAM records. Thus in
/// addition to normal fields like bases, qualities, mapping coordinates, etc.,
/// tags are used extensively to annotate records with additional
/// PacBio-specific data.
///
/// Mapping and clipping APIs are provided as well to ensure that such
/// operations "trickle down" to all data fields properly.
///
/// \sa https://samtools.github.io/hts-specs/SAMv1.pdf
///     for more information on standard %BAM data, and
///     https://github.com/PacificBiosciences/PacBioFileFormats/blob/3.0/BAM.rst
///     for more information on %PacBio %BAM fields.
///
class PBBAM_EXPORT BamRecord
{
public:
    /// \name Constructors & Related Methods
    /// \{

    BamRecord();
    BamRecord(BamHeader header);
    BamRecord(BamRecordImpl impl);
    BamRecord(const BamRecord& other);
    BamRecord(BamRecord&& other) noexcept;
    BamRecord& operator=(const BamRecord& other);
    BamRecord& operator=(BamRecord&& other) noexcept;
    virtual ~BamRecord();

    /// \}

public:
    /// \name General Data
    /// \{

    /// \returns this record's full name
    /// \sa BamRecordImpl::Name
    ///
    std::string FullName() const;

    /// \returns shared pointer to this record's associated BamHeader
    BamHeader Header() const;

    /// \returns ZMW hole number
    /// \throws if missing zm tag & record name does not contain hole number
    ///
    std::int32_t HoleNumber() const;

    /// \returns this record's LocalContextFlags
    Data::LocalContextFlags LocalContextFlags() const;

    /// \returns this record's movie name
    std::string MovieName() const;

    /// \returns "number of complete passes of the insert"
    std::int32_t NumPasses() const;

    /// \returns the record's query end position, or Sequence().length() if not
    ///          stored
    /// \note QueryEnd is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Data::Position QueryEnd() const;

    /// \returns the number of frames from start of movie to the last base of read
    ///
    std::int32_t QueryEndFrameNumber() const;

    /// \returns the record's query start position, or 0 if not stored
    ///
    /// \note QueryStart is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Data::Position QueryStart() const;

    /// \returns the number of frames from start of movie to the first base of read
    ///
    std::int32_t QueryStartFrameNumber() const;

    /// \returns this record's expected read accuracy [0, 1000]
    Data::Accuracy ReadAccuracy() const;

    /// \returns ReadGroupInfo object for this record
    ReadGroupInfo ReadGroup() const;

    /// \returns string ID of this record's read group
    ///
    /// This method should be perferred over ReadGroupBaseId() in most cases,
    /// e.g. mapping between header info.
    ///
    /// For "ID:12345678":
    ///     b.ReadGroupId()     -> "12345678"
    ///     b.ReadGroupBaseId() -> "12345678"
    ///
    /// For "ID:12345678/0--0":
    ///     b.ReadGroupId()   -> "12345678/0--0";
    ///     b.ReadGroupBaseId -> "12345678"
    ///
    /// \sa BamRecord::ReadGroupBaseId
    /// \sa ReadGroupInfo::Id
    /// \sa ReadGroupInfo::BaseId
    ///
    std::string ReadGroupId() const;

    /// \returns string base ID (stripped of optional barcode labels)
    ///
    /// ReadGroupId() should be preferred over this method in most cases. This
    /// is intended for use with hash-string or integers directly.
    ///
    /// For "ID:12345678":
    ///     b.ReadGroupId()     -> "12345678"
    ///     b.ReadGroupBaseId() -> "12345678"
    ///
    /// For "ID:12345678/0--0":
    ///     b.ReadGroupId()   -> "12345678/0--0";
    ///     b.ReadGroupBaseId -> "12345678"
    ///
    /// \sa BamRecord::ReadGroupId
    /// \sa ReadGroupInfo::Id
    /// \sa ReadGroupInfo::BaseId
    ///
    std::string ReadGroupBaseId() const;

    /// \returns integer value for this record's read group ID
    std::int32_t ReadGroupNumericId() const;

    /// \returns this scrap record's scrap region type
    VirtualRegionType ScrapRegionType() const;

    /// \returns this scrap record's scrap ZMW type
    ZmwType ScrapZmwType() const;

    /// \returns this record's average signal-to-noise for each of A, C, G,
    ///          and T
    ///
    std::vector<float> SignalToNoise() const;

    /// \returns this record's type
    /// \sa RecordType
    RecordType Type() const;

    /// \}

public:
    /// \name Mapping Data
    /// \{

    /// \returns the record's aligned end position
    ///
    /// \note AlignedEnd is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Data::Position AlignedEnd() const;

    /// \returns the record's aligned start position
    ///
    /// \note AlignedStart is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Data::Position AlignedStart() const;

    /// \returns the record's strand as a Strand enum value
    Data::Strand AlignedStrand() const;

    /// \returns the record's CIGAR data as a Cigar object
    ///
    /// \param[in] exciseAllClips   if true, remove all clipping operations
    ///                             (hard & soft) [default:false]
    ///
    Data::Cigar CigarData(bool exciseAllClips = false) const;

    /// \returns true if this record was mapped by aligner
    bool IsMapped() const;

    /// \returns this record's mapping quality. A value of 255 indicates
    ///          "unknown"
    ///
    std::uint8_t MapQuality() const;

    /// \returns the number of deleted bases (relative to reference)
    std::size_t NumDeletedBases() const;

    /// \returns the number of deletion operations (e.g. 'D' in CIGAR)
    std::size_t NumDeletionOperations() const;

    /// \returns a tuple containing NumInsertedBases (first) and NumDeletedBases
    ///         (second)
    ///
    std::pair<std::size_t, std::size_t> NumInsertedAndDeletedBases() const;

    /// \returns the number of inserted bases (relative to reference)
    std::size_t NumInsertedBases() const;

    /// \returns a tuple containing NumInsertionOperations (first) and
    ///          NumDeletionOperations (second)
    ///
    std::pair<std::size_t, std::size_t> NumInsertionAndDeletionOperations() const;

    /// \returns the number of insertion operations (e.g. 'I' in CIGAR)
    std::size_t NumInsertionOperations() const;

    /// \returns the number of matching bases (sum of '=' CIGAR op lengths)
    std::size_t NumMatches() const;

    /// \returns a tuple containing NumMatches (first) and NumMismatches
    ///         (second)
    ///
    std::pair<std::size_t, std::size_t> NumMatchesAndMismatches() const;

    /// \returns the number of mismatching bases (sum of 'X' CIGAR op lengths)
    std::size_t NumMismatches() const;

    /// \returns this record's reference ID, or -1 if unmapped.
    ///
    /// \note This is only a valid identifier within this %BAM file
    ///
    std::int32_t ReferenceId() const;

    /// \returns this record's reference name.
    ///
    /// \throws an exception if unmapped record.
    ///
    std::string ReferenceName() const;

    /// \returns the record's reference end position, or UNMAPPED_POSITION if
    ///          unmapped
    ///
    /// \note ReferenceEnd is in reference coordinates, NOT polymerase read
    ///       coordinates.
    ///
    Data::Position ReferenceEnd() const;

    /// \returns the record's reference start position, or UNMAPPED_POSITION if
    ///          unmapped
    ///
    /// \note ReferenceStart is in reference coordinates, NOT polymerase read
    ///       coordinates.
    ///
    Data::Position ReferenceStart() const;

    /// \}

public:
    /// \name Barcode Data
    /// \{

    /// \returns forward barcode id
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    std::int16_t BarcodeForward() const;

    /// \returns barcode call confidence (Phred-scaled posterior probability
    ///          of correct barcode call)
    ///
    /// \sa HasBarcodeQuality
    ///
    std::uint8_t BarcodeQuality() const;

    /// \returns reverse barcode id
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    std::int16_t BarcodeReverse() const;

    /// \returns the forward and reverse barcode ids
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    std::pair<std::int16_t, std::int16_t> Barcodes() const;

    /// \}

public:
    /// \name Segment Read Data

    /// \returns true if segmented read
    ///
    /// \note Queries the record's read group, not the tags.
    ///
    bool IsSegment() const;

    /// \returns segment read index
    ///
    /// \throws std::runtime_error if segment read data is absent or malformed.
    /// \sa HasSegmentIndex
    ///
    std::int32_t SegmentIndex() const;

    /// \returns index of left adapater
    ///
    /// \throws std::runtime_error if segment read data is absent or malformed.
    /// \sa HasSegmentLeadingAdapterIndex
    ///
    std::int32_t SegmentLeftAdapterIndex() const;

    /// \returns index of right adapater
    ///
    /// \throws std::runtime_error if segment read data is absent or malformed.
    /// \sa HasSegmentTrailingAdapterIndex
    ///
    std::int32_t SegmentRightAdapterIndex() const;

    /// \returns segment read supplemental data, decoded to JSON
    ///
    /// \throws std::runtime_error if segment read data is absent or malformed.
    /// \sa HasSegmentSupplementalData
    ///
    JSON::Json SegmentSupplementalData() const;

    /// \}

public:
    /// \name Auxiliary Data Queries
    /// \{

    /// \returns true if this record has AltLabelQV data
    bool HasAltLabelQV() const;

    /// \returns true if this record has AltLabelTag data
    bool HasAltLabelTag() const;

    /// \returns true if this record has Barcode data
    bool HasBarcodes() const;

    /// \returns true is this record has BarcodeQuality data
    bool HasBarcodeQuality() const;

    /// \returns true if this record has DeletionQV data
    bool HasDeletionQV() const;

    /// \returns true if this record has DeletionTag data
    bool HasDeletionTag() const;

    /// \returns true if this record has forward IPD data
    bool HasForwardIPD() const;

    /// \returns true if this record has forward pulse width data
    bool HasForwardPulseWidth() const;

    /// \returns true if this record has a HoleNumber
    bool HasHoleNumber() const;

    /// \returns true if this record has InsertionQV data
    bool HasInsertionQV() const;

    /// \returns true if this record has IPD data
    bool HasIPD() const;

    /// \returns true if this record has LabelQV data
    bool HasLabelQV() const;

    /// \returns true if this record has LocalContextFlags (absent in CCS)
    bool HasLocalContextFlags() const;

    /// \returns true if this record has MergeQV data
    bool HasMergeQV() const;

    /// \returns true if this record has NumPasses data
    bool HasNumPasses() const;

    /// \returns true if this record has Pkmean data
    bool HasPkmean() const;

    /// \returns true if this record has Pkmid data
    bool HasPkmid() const;

    /// \returns true if this record has Pkmean2 data
    bool HasPkmean2() const;

    /// \returns true if this record has Pkmid2 data
    bool HasPkmid2() const;

    /// \returns true if this record has PreBaseFrames aka IPD data
    bool HasPreBaseFrames() const;

    /// \returns true if this record has PrePulseFrames data
    bool HasPrePulseFrames() const;

    /// \returns true if this record has PulseCall data
    bool HasPulseCall() const;

    /// \returns true if this record has PulseCallWidth data
    bool HasPulseCallWidth() const;

    /// \returns true if this record has PulseExclusion data
    bool HasPulseExclusion() const;

    /// \returns true if this record has PulseMergeQV data
    bool HasPulseMergeQV() const;

    /// \returns true if this record has PulseWidth data
    bool HasPulseWidth() const;

    /// \returns true if this record has ReadAccuracyTag data
    bool HasReadAccuracy() const;

    /// \returns true if this record has QueryEnd data
    bool HasQueryEnd() const;

    /// \returns true if this record has QueryEnd data
    bool HasQueryEndFrameNumber() const;

    /// \returns true if this record has QueryStart data
    bool HasQueryStart() const;

    /// \returns true if this record has QueryStartFrameNumber data
    bool HasQueryStartFrameNumber() const;

    /// \returns true if this record has reverse IPD data
    bool HasReverseIPD() const;

    /// \returns true if this record has reverse pulse width data
    bool HasReversePulseWidth() const;

    /// \returns true if this record has ScrapRegionType data (only in SCRAP)
    bool HasScrapRegionType() const;

    /// \returns true if this record has scrap ZMW type data (only in SCRAP)
    bool HasScrapZmwType() const;

    /// \returns true if this record has segment index
    bool HasSegmentIndex() const;

    /// \returns true if this record has segment's left adapter index
    bool HasSegmentLeftAdapterIndex() const;

    /// \returns true if this record has segment's right adapter index
    bool HasSegmentRightAdapterIndex() const;

    /// \returns true if this record has segment supplemental data
    bool HasSegmentSupplementalData() const;

    /// \returns true if this record has signal-to-noise data (absent in
    ///          POLYMERASE)
    ///
    bool HasSignalToNoise() const;

    /// \returns true if this record has StartFrame data
    bool HasStartFrame() const;

    /// \returns true if this record has SubstitutionQV data
    bool HasSubstitutionQV() const;

    /// \returns true if this record has SubstitutionTag data
    bool HasSubstitutionTag() const;

    /// \}

public:
    /// \name Sequence & Tag Data
    /// \{

    /// \brief Fetches this record's AltLabelTag values ("pt" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new gap chars will be '-' and padding chars will be '*'.
    ///
    /// \param[in] orientation      Orientation of output.
    ///
    /// \returns AltLabelTags string
    ///
    std::string AltLabelTag(Data::Orientation orientation = Data::Orientation::NATIVE,
                            bool aligned = false, bool exciseSoftClips = false,
                            PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's DeletionTag values ("dt" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new gap chars will be '-' and padding chars will be '*'.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns DeletionTag string
    ///
    std::string DeletionTag(Data::Orientation orientation = Data::Orientation::NATIVE,
                            bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's DNA sequence (SEQ field).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new gap chars will be '-' and padding chars will be '*'.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns sequence string
    ///
    std::string Sequence(Data::Orientation orientation = Data::Orientation::NATIVE,
                         bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's SubstitutionTag values ("st" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new gap chars will be '-' and padding chars will be '*'.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns SubstitutionTags string
    ///
    std::string SubstitutionTag(Data::Orientation orientation = Data::Orientation::NATIVE,
                                bool aligned = false, bool exciseSoftClips = false) const;

    /// \}

public:
    /// \name Quality Data
    /// \{

    /// \brief Fetches this record's AltLabelQV values ("pv" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    ///
    /// \returns AltLabelQV as QualityValues object
    ///
    Data::QualityValues AltLabelQV(Data::Orientation orientation = Data::Orientation::NATIVE,
                                   bool aligned = false, bool exciseSoftClips = false,
                                   PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's DeletionQV values ("dq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns DeletionQV as QualityValues object
    ///
    Data::QualityValues DeletionQV(Data::Orientation orientation = Data::Orientation::NATIVE,
                                   bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's InsertionQV values ("iq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns InsertionQVs as QualityValues object
    ///
    Data::QualityValues InsertionQV(Data::Orientation orientation = Data::Orientation::NATIVE,
                                    bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's LabelQV values ("pq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    ///
    /// \returns LabelQV as QualityValues object
    ///
    Data::QualityValues LabelQV(Data::Orientation orientation = Data::Orientation::NATIVE,
                                bool aligned = false, bool exciseSoftClips = false,
                                PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's MergeQV values ("mq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns MergeQV as QualityValues object
    ///
    Data::QualityValues MergeQV(Data::Orientation orientation = Data::Orientation::NATIVE,
                                bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches  this record's %BAM quality values (QUAL field).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns %BAM qualities as QualityValues object
    ///
    Data::QualityValues Qualities(Data::Orientation orientation = Data::Orientation::NATIVE,
                                  bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's SubstitutionQV values ("sq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns SubstitutionQV as QualityValues object
    ///
    Data::QualityValues SubstitutionQV(Data::Orientation orientation = Data::Orientation::NATIVE,
                                       bool aligned = false, bool exciseSoftClips = false) const;

    /// \}

public:
    /// \name Pulse Data
    /// \{

    /// \brief Fetches this record's forward IPD values ("fi" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new frames will have a value of 0;
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns forward IPD as Frames object
    ///
    Data::Frames ForwardIPD(Data::Orientation orientation = Data::Orientation::NATIVE,
                            bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's forward pulse width values ("fp" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new frames will have a value of 0;
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns foward PW as Frames object
    ///
    Data::Frames ForwardPulseWidth(Data::Orientation orientation = Data::Orientation::NATIVE,
                                   bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's IPD values ("ip" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new frames will have a value of 0;
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns IPD as Frames object
    ///
    Data::Frames IPD(Data::Orientation orientation = Data::Orientation::NATIVE,
                     bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's IPD values ("ip" tag), but does not upscale.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns IPD as Frames object
    ///
    Data::Frames IPDRaw(Data::Orientation orientation = Data::Orientation::NATIVE) const;

    /// \brief Fetches this record's Pkmean values ("pa" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmean as vector<float> object
    ///
    std::vector<float> Pkmean(Data::Orientation orientation = Data::Orientation::NATIVE,
                              bool aligned = false, bool exciseSoftClips = false,
                              PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's Pkmid values ("pm" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmid as vector<float> object
    ///
    std::vector<float> Pkmid(Data::Orientation orientation = Data::Orientation::NATIVE,
                             bool aligned = false, bool exciseSoftClips = false,
                             PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's Pkmean2 values ("pi" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmean as vector<float> object
    ///
    std::vector<float> Pkmean2(Data::Orientation orientation = Data::Orientation::NATIVE,
                               bool aligned = false, bool exciseSoftClips = false,
                               PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's Pkmid2 values ("ps" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmid as vector<float> object
    ///
    std::vector<float> Pkmid2(Data::Orientation orientation = Data::Orientation::NATIVE,
                              bool aligned = false, bool exciseSoftClips = false,
                              PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PreBaseFrames aka IPD values ("ip" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new frames will have a value of 0;
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns IPD as Frames object
    ///
    Data::Frames PreBaseFrames(Data::Orientation orientation = Data::Orientation::NATIVE,
                               bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's PrePulseFrames values ("pd" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PrePulseFrames as Frames object
    ///
    Data::Frames PrePulseFrames(Data::Orientation orientation = Data::Orientation::NATIVE,
                                bool aligned = false, bool exciseSoftClips = false,
                                PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PulseCall values ("pc" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseCalls string
    ///
    std::string PulseCall(Data::Orientation orientation = Data::Orientation::NATIVE,
                          bool aligned = false, bool exciseSoftClips = false,
                          PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PulseCallWidth values ("px" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseCallWidth as Frames object
    ///
    Data::Frames PulseCallWidth(Data::Orientation orientation = Data::Orientation::NATIVE,
                                bool aligned = false, bool exciseSoftClips = false,
                                PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PulseExclusionReason values ("pe" tag).
    ///
    /// \returns vector of pulse exclusion reason value
    ///
    std::vector<BAM::PulseExclusionReason> PulseExclusionReason(
        Data::Orientation orientation = Data::Orientation::NATIVE, bool aligned = false,
        bool exciseSoftClips = false, PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetch this record's PulseMergeQV values ("pg" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseMergeQV as QualityValues object
    ///
    Data::QualityValues PulseMergeQV(Data::Orientation orientation = Data::Orientation::NATIVE,
                                     bool aligned = false, bool exciseSoftClips = false,
                                     PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PulseWidth values ("pw" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new frames will have a value of 0.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns PulseWidths as Frames object
    ///
    Data::Frames PulseWidth(Data::Orientation orientation = Data::Orientation::NATIVE,
                            bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's PulseWidth values ("pw" tag), but does not
    ///        upscale.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseWidth as Frames object
    ///
    Data::Frames PulseWidthRaw(Data::Orientation orientation = Data::Orientation::NATIVE,
                               bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's reverse IPD values ("ri" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new frames will have a value of 0;
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns reverse IPD as Frames object
    ///
    Data::Frames ReverseIPD(Data::Orientation orientation = Data::Orientation::NATIVE,
                            bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's reverse pulse width values ("rp" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new frames will have a value of 0.
    ///
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns reverse PulseWidths as Frames object
    ///
    Data::Frames ReversePulseWidth(Data::Orientation orientation = Data::Orientation::NATIVE,
                                   bool aligned = false, bool exciseSoftClips = false) const;

    /// \brief Fetches this record's StartFrame values ("sf" tag).
    ///
    /// \param[in] orientation     Orientation of output
    ///
    /// \returns StartFrame as std::uint32_t vector
    ///
    std::vector<std::uint32_t> StartFrame(Data::Orientation orientation = Data::Orientation::NATIVE,
                                          bool aligned = false, bool exciseSoftClips = false,
                                          PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \}

public:
    /// \name Low-Level Access & Operations
    /// \{

    /// \warning This method should be considered temporary and avoided as much
    ///          as possible. Direct access to the internal object is likely to
    ///          disappear as BamRecord interface matures.
    ///
    /// \returns const reference to underlying BamRecordImpl object
    ///
    const BamRecordImpl& Impl() const;

    /// \warning This method should be considered temporary and avoided as much
    ///          as possible. Direct access to the internal object is likely to
    ///          disappear as BamRecord interface matures.
    ///
    /// \returns reference to underlying BamRecordImpl object
    ///
    BamRecordImpl& Impl();

    /// \}

public:
    /// \name General Data
    /// \{

    /// \brief Sets this record's ZMW hole number.
    ///
    /// \param[in] holeNumber
    /// \returns reference to this record
    ///
    BamRecord& HoleNumber(std::int32_t holeNumber);

    /// \brief Sets this record's local context flags
    ///
    /// \param[in] flags
    /// \returns reference to this record
    ///
    BamRecord& LocalContextFlags(Data::LocalContextFlags flags);

    /// \brief Sets this record's "number of complete passes of the insert".
    ///
    /// \param[in] numPasses
    /// \returns reference to this record
    ///
    BamRecord& NumPasses(std::int32_t numPasses);

    /// \brief Sets this record's query end position.
    ///
    /// \note Changing this will modify the name of non-CCS records.
    ///
    /// \param[in] pos
    /// \returns reference to this record
    ///
    BamRecord& QueryEnd(Data::Position pos);

    /// \brief Sets this record's query end frame number
    ///
    /// \param[in] frame number
    /// \returns reference to this record
    ///
    BamRecord& QueryEndFrameNumber(std::int32_t frameNumber);

    /// \brief Sets this record's query start position.
    ///
    /// \note Changing this will modify the name of non-CCS records.
    ///
    /// \param[in] pos
    /// \returns reference to this record
    ///
    BamRecord& QueryStart(Data::Position pos);

    /// \brief Sets this record's query start frame number
    ///
    /// \param[in] frame number
    /// \returns reference to this record
    ///
    BamRecord& QueryStartFrameNumber(std::int32_t frameNumber);

    /// \brief Sets this record's expected read accuracy [0, 1000]
    ///
    /// \param[in] accuracy
    /// \returns reference to this record
    ///
    BamRecord& ReadAccuracy(const Data::Accuracy& accuracy);

    /// \brief Attaches this record to the provided read group, changing the
    ///        record name & 'RG' tag.
    ///
    /// \param[in] rg
    /// \returns reference to this record
    ///
    BamRecord& ReadGroup(const ReadGroupInfo& rg);

    /// \brief Attaches this record to the provided read group, changing the
    ///        record name & 'RG' tag.
    ///
    /// \param[in] id
    /// \returns reference to this record
    ///
    BamRecord& ReadGroupId(const std::string& id);

    /// \brief Sets this scrap record's ScrapRegionType
    ///
    /// \param[in] type
    /// \returns reference to this record
    ///
    BamRecord& ScrapRegionType(VirtualRegionType type);

    /// \brief Sets this scrap record's ScrapRegionType
    ///
    /// \param[in] type character equivalent of VirtualRegionType
    /// \returns reference to this record
    ///
    BamRecord& ScrapRegionType(char type);

    /// \brief Sets this scrap record's ScrapZmwType
    ///
    /// \param[in] type
    /// \returns reference to this record
    ///
    BamRecord& ScrapZmwType(ZmwType type);

    /// \brief Sets this scrap record's ScrapZmwType
    ///
    /// \param[in] type character equivalent of ZmwType
    /// \returns reference to this record
    ///
    BamRecord& ScrapZmwType(char type);

    /// \brief Sets this record's average signal-to-noise in each of A, C, G,
    ///        and T
    ///
    /// \param[in] snr average signal-to-noise of A, C, G, and T (in this order)
    /// \returns reference to this record
    ///
    BamRecord& SignalToNoise(const std::vector<float>& snr);

    /// \}

public:
    /// \name Barcode Data
    /// \{

    /// \brief Sets this record's barcode IDs ('bc' tag)
    ///
    /// \param[in] barcodeIds
    /// \returns reference to this record
    ///
    BamRecord& Barcodes(const std::pair<std::int16_t, std::int16_t>& barcodeIds);

    /// \brief Sets this record's barcode quality ('bq' tag)
    ///
    /// \param[in] quality Phred-scaled confidence call
    /// \returns reference to this record
    ///
    BamRecord& BarcodeQuality(std::uint8_t quality);

    /// \}

public:
    /// \name Segment Data

    /// \brief Sets this record's segment index ('di' tag)
    ///
    /// \param[in] index    0-based index of this segment within its source read
    /// \returns reference to this record
    ///
    BamRecord& SegmentIndex(std::int32_t index);

    /// \returns Sets this segment's left adapter index ('dl' tag)
    ///
    /// \param[in] index    0-based index of segment left adapter
    /// \returns reference to this record
    ///
    BamRecord& SegmentLeftAdapterIndex(std::int32_t index);

    /// \returns Sets this segment's right adapter index ('dr' tag)
    ///
    /// \param[in] index    0-based index of segment right adapter
    /// \returns reference to this record
    ///
    BamRecord& SegmentRightAdapterIndex(std::int32_t index);

    /// \returns Sets this segment's supplemental data ('ds' tag)
    ///
    /// \param[in] data     JSON object (not encoded)
    /// \returns reference to this record
    ///
    BamRecord& SegmentSupplementalData(const JSON::Json& data);

    /// \}

public:
    /// \name Sequence & Tag Data
    /// \{

    /// \brief Sets this record's AltLabelTag values ("at" tag).
    ///
    /// \param[in] tags
    /// \returns reference to this record
    ///
    BamRecord& AltLabelTag(const std::string& tags);

    /// \brief Sets this record's DeletionTag values ("dt" tag).
    ///
    /// \param[in] tags
    /// \returns reference to this record
    ///
    BamRecord& DeletionTag(const std::string& tags);

    /// \brief Sets this record's SubstitutionTag values ("st" tag).
    ///
    /// \param[in] tags
    /// \returns reference to this record
    ///
    BamRecord& SubstitutionTag(const std::string& tags);

    /// \}

public:
    /// \name Quality Data
    /// \{

    /// \brief Sets this record's AltLabelQV values ("pv" tag).
    ///
    /// \param[in] altLabelQVs
    /// \returns reference to this record
    ///
    BamRecord& AltLabelQV(const Data::QualityValues& altLabelQVs);

    /// \brief Sets this record's DeletionQV values ("dq" tag).
    ///
    /// \param[in] deletionQVs
    /// \returns reference to this record
    ///
    BamRecord& DeletionQV(const Data::QualityValues& deletionQVs);

    /// \brief Sets this record's InsertionQV values ("iq" tag).
    ///
    /// \param[in] insertionQVs
    /// \returns reference to this record
    ///
    BamRecord& InsertionQV(const Data::QualityValues& insertionQVs);

    /// \brief Sets this record's LabelQV values ("pq" tag).
    ///
    /// \param[in] labelQVs
    /// \returns reference to this record
    ///
    BamRecord& LabelQV(const Data::QualityValues& labelQVs);

    /// \brief Sets this record's MergeQV values ("mq" tag).
    ///
    /// \param[in] mergeQVs
    /// \returns reference to this record
    ///
    BamRecord& MergeQV(const Data::QualityValues& mergeQVs);

    /// \brief Sets this record's SubstitutionQV values ("sq" tag).
    ///
    /// \param[in] substitutionQVs
    /// \returns reference to this record
    ///
    BamRecord& SubstitutionQV(const Data::QualityValues& substitutionQVs);

    /// \}

public:
    /// \name Pulse Data
    /// \{

    /// \brief Sets this record's forward IPD values ("fi" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& ForwardIPD(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's forward pulse width values ("fp" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& ForwardPulseWidth(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's IPD values ("ip" tag).
    ///
    /// \deprecated since v1.3.0. Use the FrameCodec overload instead
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    PBBAM_DEPRECATED_FRAMES BamRecord& IPD(const Data::Frames& frames, FrameEncodingType encoding);

    /// \brief Sets this record's IPD values ("ip" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& IPD(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's Pkmean values ("pm" tag).
    ///
    /// \param[in] photons
    /// \returns reference to this record
    ///
    BamRecord& Pkmean(const std::vector<float>& photons);

    /// \brief Sets this record's Pkmean values ("pm" tag).
    ///
    /// \param[in] encodedPhotons
    /// \returns reference to this record
    ///
    BamRecord& Pkmean(const std::vector<std::uint16_t>& encodedPhotons);

    /// \brief Sets this record's Pkmid values ("pa" tag).
    ///
    /// \param[in] photons
    /// \returns reference to this record
    ///
    BamRecord& Pkmid(const std::vector<float>& photons);

    /// \brief Sets this record's Pkmid values ("pa" tag).
    ///
    /// \param[in] encodedPhotons
    /// \returns reference to this record
    ///
    BamRecord& Pkmid(const std::vector<std::uint16_t>& encodedPhotons);

    /// \brief Sets this record's Pkmean2 values ("ps" tag).
    ///
    /// \param[in] photons
    /// \returns reference to this record
    ///
    BamRecord& Pkmean2(const std::vector<float>& photons);

    /// \brief Sets this record's Pkmean2 values ("ps" tag).
    ///
    /// \param[in] encodedPhotons
    /// \returns reference to this record
    ///
    BamRecord& Pkmean2(const std::vector<std::uint16_t>& encodedPhotons);

    /// \brief Sets this record's Pkmid2 values ("pi" tag).
    ///
    /// \param[in] photons
    /// \returns reference to this record
    ///
    BamRecord& Pkmid2(const std::vector<float>& photons);

    /// \brief Sets this record's Pkmid2 values ("pi" tag).
    ///
    /// \param[in] encodedPhotons
    /// \returns reference to this record
    ///
    BamRecord& Pkmid2(const std::vector<std::uint16_t>& encodedPhotons);

    /// \brief Sets this record's PreBaseFrames aka IPD values ("ip" tag).
    ///
    /// \deprecated since v1.3.0. Use the FrameCodec overload instead
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    PBBAM_DEPRECATED_FRAMES BamRecord& PreBaseFrames(const Data::Frames& frames,
                                                     FrameEncodingType encoding);

    /// \brief Sets this record's PreBaseFrames aka IPD values ("ip" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PreBaseFrames(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's PrePulseFrames values ("pd" tag).
    ///
    /// \deprecated since v1.3.0. Use the FrameCodec overload instead
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    PBBAM_DEPRECATED_FRAMES BamRecord& PrePulseFrames(const Data::Frames& frames,
                                                      FrameEncodingType encoding);

    /// \brief Sets this record's PrePulseFrames values ("pd" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PrePulseFrames(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's PulseCall values ("pc" tag).
    ///
    /// \param[in] tags
    /// \returns reference to this record
    ///
    BamRecord& PulseCall(const std::string& tags);

    /// \brief Sets this record's PulseCallWidth values ("px" tag).
    ///
    /// \deprecated since v1.3.0. Use the FrameCodec overload instead
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    PBBAM_DEPRECATED_FRAMES BamRecord& PulseCallWidth(const Data::Frames& frames,
                                                      FrameEncodingType encoding);

    /// \brief Sets this record's PulseCallWidth values ("px" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PulseCallWidth(const Data::Frames& frames, Data::FrameCodec encoding);

    ///
    /// \\brief Sets this record's PulseExclusionReason values ("pe" tag).
    /// \param[in] reasons
    /// \return reference to this record
    ///
    BamRecord& PulseExclusionReason(const std::vector<BAM::PulseExclusionReason>& reasons);

    /// \brief Sets this record's PulseMergeQV values ("pg" tag).
    ///
    /// \param[in] pulseMergeQVs
    /// \returns reference to this record
    ///
    BamRecord& PulseMergeQV(const Data::QualityValues& pulseMergeQVs);

    /// \brief Sets this record's PulseWidth values ("pw" tag).
    ///
    /// \deprecated since v1.3.0. Use the FrameCodec overload instead
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    PBBAM_DEPRECATED_FRAMES BamRecord& PulseWidth(const Data::Frames& frames,
                                                  FrameEncodingType encoding);

    /// \brief Sets this record's PulseWidth values ("pw" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PulseWidth(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's reverse IPD values ("ri" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& ReverseIPD(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's reverse pulse width values ("rp" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& ReversePulseWidth(const Data::Frames& frames, Data::FrameCodec encoding);

    /// \brief Sets this record's StartFrame values ("sf" tag).
    ///
    /// \param[in] startFrame
    /// \returns reference to this record
    ///
    BamRecord& StartFrame(const std::vector<std::uint32_t>& startFrame);

    /// \}

public:
    /// \name Low-Level Access & Operations
    /// \{

    /// \brief Resets cached aligned start/end.
    ///
    /// \note This method should not be needed in most client code. It exists
    ///       primarily as a hook for internal reading loops (queries, index
    ///       build, etc.) It's essentially a workaround and will likely be
    ///       removed from the API.
    ///
    void ResetCachedPositions() const;

    /// \brief Resets cached aligned start/end.
    ///
    /// \note This method should not be needed in most client code. It exists
    ///       primarily as a hook for internal reading loops (queries, index
    ///       build, etc.) It's essentially a workaround and will likely be
    ///       removed from the API.
    ///
    void ResetCachedPositions();

    /// \brief Updates the record's name (BamRecord::FullName) to reflect
    ///        modifications to name components (movie name, ZMW hole number,
    ///        etc.)
    ///
    void UpdateName();

    /// \}

public:
    /// \name Pulse Data
    /// \{

    static const float photonFactor;

    static std::vector<std::uint16_t> EncodePhotons(const std::vector<float>& data);

    /// \}

public:
    /// \name [Mapped]Read conversion
    /// \{

    ///
    /// \return Data::Read representation of this record
    ///
    Data::Read ToRead(std::string model = "") const;

    ///
    /// \return Data::MappedRead representation of this record
    ///
    /// \throws std::runtime_error if record is unmapped
    ///
    Data::MappedRead ToMappedRead(std::string model = "", Data::Position startOffset = 0,
                                  bool pinStart = false, bool pinEnd = false) const;

    /// \}

public:
    /// \name Clipping & Mapping
    /// \{

    /// Creates a copied record from input, with clipping applied
    static BamRecord Clipped(const BamRecord& input, ClipType clipType, Data::Position start,
                             Data::Position end, bool exciseFlankingInserts = false);

    /// Creates a copied record from input, with mapping applied
    static BamRecord Mapped(const BamRecord& input, std::int32_t referenceId,
                            Data::Position refStart, Data::Strand strand, const Data::Cigar& cigar,
                            std::uint8_t mappingQuality);

    /// Splits the (5mC) basemods `Mm` and `Ml` tags
    struct SplitBasemods
    {
        std::vector<std::int32_t> LeadingSeparatingC;
        std::vector<std::uint8_t> LeadingQuals;

        std::vector<std::int32_t> RetainedSeparatingC;
        std::vector<std::uint8_t> RetainedQuals;

        std::vector<std::int32_t> TrailingSeparatingC;
        std::vector<std::uint8_t> TrailingQuals;

        std::int32_t PrefixLostBases{0};

        static std::vector<std::int32_t> SplitBasemodsString(const std::string& str);

        static std::string SeparatingCToString(const std::vector<std::int32_t>& vec);
    };
    static SplitBasemods ClipBasemodsTag(const std::string& seq,
                                         const std::string& oldBasemodsString,
                                         const std::vector<std::uint8_t>& basemodsQVs,
                                         std::size_t clipFrom, std::size_t clipLength);

    /// Splits subread pileup tags 'sa', 'sm' and 'sx'
    struct SplitSubreadPileup
    {
        std::vector<std::uint16_t> LeadingCoverage{};
        std::vector<std::uint8_t> LeadingMatches{};
        std::vector<std::uint8_t> LeadingMismatches{};

        std::vector<std::uint16_t> RetainedCoverage{};
        std::vector<std::uint8_t> RetainedMatches{};
        std::vector<std::uint8_t> RetainedMismatches{};

        std::vector<std::uint16_t> TrailingCoverage{};
        std::vector<std::uint8_t> TrailingMatches{};
        std::vector<std::uint8_t> TrailingMismatches{};

        std::int32_t LostPrefixBases{0};
        std::int32_t LostSuffixBases{0};
        std::int32_t LostCoverage{-1};

        bool operator==(const SplitSubreadPileup&) const = default;
    };
    static SplitSubreadPileup ClipSubreadPileupTags(
        std::size_t sequenceLength, const std::vector<std::uint16_t>& runLengthEncodedCoverage,
        const std::vector<std::uint8_t>& matches, const std::vector<std::uint8_t>& mismatches,
        std::size_t clipFrom, std::size_t clipLength);

    /// Applies clipping to this record
    BamRecord& Clip(ClipType clipType, Data::Position start, Data::Position end,
                    bool exciseFlankingInserts = false);

    /// Creates a copied record from this one, with clipping applied
    BamRecord Clipped(ClipType clipType, Data::Position start, Data::Position end,
                      bool exciseFlankingInserts = false) const;

    /// Applies mapping to this record
    BamRecord& Map(std::int32_t referenceId, Data::Position refStart, Data::Strand strand,
                   const Data::Cigar& cigar, std::uint8_t mappingQuality);

    /// Creates a copied record from this one, with mapping applied
    BamRecord Mapped(std::int32_t referenceId, Data::Position refStart, Data::Strand strand,
                     const Data::Cigar& cigar, std::uint8_t mappingQuality) const;
    /// \}

    ///
    /// \returns estimated number of bytes used by this record
    ///
    /// \warning The actual usage is heavily implementation-dependent, w.r.t.
    ///          data structure layout and alignment. A general estimate is
    ///          provided here, but no guarantee can be made.
    ///
    int EstimatedBytesUsed() const noexcept;

private:
    BamRecordImpl impl_;

public:
    /// public & mutable so that queries can directly set the header info,
    /// even on a record that is const from client code's perspective
    mutable BamHeader header_;

private:
    /// \internal
    /// cached positions (mutable to allow lazy-calc in const methods)
    mutable Data::Position alignedStart_ = Data::UNMAPPED_POSITION;
    mutable Data::Position alignedEnd_ = Data::UNMAPPED_POSITION;

private:
    /// \internal
    /// pulse to bam mapping cache
    mutable std::unique_ptr<Pulse2BaseCache> p2bCache_;

public:
    /// clips the PacBio tags to a specified length
    void ClipTags(std::size_t clipPos, std::size_t clipLength);

private:
    ///\internal
    /// clipping methods

    void ClipFields(std::size_t clipPos, std::size_t clipLength);

    BamRecord& ClipToQuery(Data::Position start, Data::Position end);
    BamRecord& ClipToReference(Data::Position start, Data::Position end,
                               bool exciseFlankingInserts);

private:
    ///\internal
    /// raw tag data fetching

    // sequence tags
    std::string FetchBasesRaw(BamRecordTag tag) const;
    std::string FetchBases(BamRecordTag tag,
                           Data::Orientation orientation = Data::Orientation::NATIVE,
                           bool aligned = false, bool exciseSoftClips = false,
                           PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // frame tags
    Data::Frames FetchFramesRaw(BamRecordTag tag) const;
    Data::Frames FetchFrames(BamRecordTag tag,
                             Data::Orientation orientation = Data::Orientation::NATIVE,
                             bool aligned = false, bool exciseSoftClips = false,
                             PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // pulse tags
    std::vector<float> FetchPhotonsRaw(BamRecordTag tag) const;
    std::vector<float> FetchPhotons(BamRecordTag tag,
                                    Data::Orientation orientation = Data::Orientation::NATIVE,
                                    bool aligned = false, bool exciseSoftClips = false,
                                    PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // QV tags
    Data::QualityValues FetchQualitiesRaw(BamRecordTag tag) const;
    Data::QualityValues FetchQualities(BamRecordTag tag,
                                       Data::Orientation orientation = Data::Orientation::NATIVE,
                                       bool aligned = false, bool exciseSoftClips = false,
                                       PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // UInt tags (e.g. start frame)
    //
    // TODO (DB): clean this up w.r.t FetchUInt8s
    //
    std::vector<std::uint32_t> FetchUInt32sRaw(BamRecordTag tag) const;
    std::vector<std::uint32_t> FetchUInt32s(
        BamRecordTag tag, Data::Orientation orientation = Data::Orientation::NATIVE,
        bool aligned = false, bool exciseSoftClips = false,
        PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // UInt tags (e.g. pulse exclusion)
    //
    // ODO (DB): clean this up w.r.t FetchUInt32s
    //
    std::vector<std::uint8_t> FetchUInt8sRaw(BamRecordTag tag) const;
    std::vector<std::uint8_t> FetchUInt8s(BamRecordTag tag,
                                          Data::Orientation orientation = Data::Orientation::NATIVE,
                                          bool aligned = false, bool exciseSoftClips = false,
                                          PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

private:
    ///\internal
    /// marked const to allow calling from const methods
    /// but updates our mutable cached values
    void CalculateAlignedPositions() const;
    void CalculatePulse2BaseCache() const;

    friend class BamRecordMemory;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BAMRECORD_H
