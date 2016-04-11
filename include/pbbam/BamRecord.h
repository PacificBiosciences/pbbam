// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// File Description
/// \file BamRecord.h
/// \brief Defines the BamRecord & BamRecordView classes.
//
// Author: Derek Barnett

#ifndef BAMRECORD_H
#define BAMRECORD_H

#include "pbbam/Accuracy.h"
#include "pbbam/Frames.h"
#include "pbbam/BamRecordImpl.h"
#include "pbbam/BamHeader.h"
#include "pbbam/LocalContextFlags.h"
#include "pbbam/Orientation.h"
#include "pbbam/ReadGroupInfo.h"
#include "pbbam/Strand.h"
#include "pbbam/QualityValues.h"
#include "pbbam/virtual/VirtualRegionType.h"
#include "pbbam/ZmwType.h"
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal { class BamRecordMemory; }

/// \brief This enum defines the modes supported by BamRecord clipping
///        operations.
///
/// Methods like BamRecord::Clip accept Position parameters - which may be in
/// either polymerase or reference coorindates. Using this enum as a flag
/// indicates how the positions should be interpreted.
///
enum class ClipType
{
    CLIP_NONE           ///< No clipping will be performed.
  , CLIP_TO_QUERY       ///< Clipping positions are in polymerase coordinates.
  , CLIP_TO_REFERENCE   ///< Clipping positions are in genomic coordinates.
};

/// \brief This enum defines the possible PacBio BAM record types.
///
/// \sa ReadGroupInfo::ReadType
///
enum class RecordType
{
    ZMW         ///< Polymerase read
  , HQREGION    ///< High-quality region
  , SUBREAD     ///< Subread (
  , CCS         ///< Circular consensus sequence
  , SCRAP       ///< Additional sequence (barcodes, adapters, etc.)
  , UNKNOWN     ///< Unknown read type

  , POLYMERASE = ZMW ///< \deprecated as of PacBio BAM spec v 3.0.4 (use RecordType::ZMW instead)
};

/// \brief This enum defines the possible encoding modes used in Frames data
/// (e.g. BamRecord::IPD or BamRecord::PulseWidth).
///
/// The LOSSY mode is the default in production output; LOSSLESS mode
/// being used primarily for internal applications.
///
/// \sa https://github.com/PacificBiosciences/PacBioFileFormats/blob/3.0/BAM.rst
///     for more information on pulse frame encoding schemes.
///
enum class FrameEncodingType
{
    LOSSY       ///< 8-bit compression (using CodecV1) of frame data
  , LOSSLESS    ///< 16-bit native frame data
};

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

    BamRecord(void);
    BamRecord(const BamHeader& header);
    BamRecord(const BamRecordImpl& impl);
    BamRecord(BamRecordImpl&& impl);
    BamRecord(const BamRecord& other);
    BamRecord(BamRecord&& other);
    BamRecord& operator=(const BamRecord& other);
    BamRecord& operator=(BamRecord&& other);
    virtual ~BamRecord(void);

    /// \}

public:
    /// \name General Data
    /// \{

    /// \returns this record's full name
    /// \sa BamRecordImpl::Name
    ///
    std::string FullName(void) const;

    /// \returns shared pointer to this record's associated BamHeader
    BamHeader Header(void) const;

    /// \returns ZMW hole number
    /// \throws if missing zm tag & record name does not contain hole number
    ///
    int32_t HoleNumber(void) const;

    /// \returns this record's LocalContextFlags
    PacBio::BAM::LocalContextFlags LocalContextFlags(void) const;

    /// \returns this record's movie name
    std::string MovieName(void) const;

    /// \returns "number of complete passes of the insert"
    int32_t NumPasses(void) const;

    /// \returns the record's query end position, or Sequence().length() if not
    ///          stored
    /// \note QueryEnd is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Position QueryEnd(void) const;

    /// \returns the record's query start position, or 0 if not stored
    ///
    /// \note QueryStart is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Position QueryStart(void) const;

    /// \returns this record's expected read accuracy [0, 1000]
    Accuracy ReadAccuracy(void) const;

    /// \returns ReadGroupInfo object for this record
    ReadGroupInfo ReadGroup(void) const;

    /// \returns string ID of this record's read group
    /// \sa ReadGroupInfo::Id
    ///
    std::string ReadGroupId(void) const;

    /// \returns integer value for this record's read group ID
    int32_t ReadGroupNumericId(void) const;

    /// \returns this scrap record's scrap region type
    VirtualRegionType ScrapRegionType(void) const;

    /// \returns this scrap record's scrap ZMW type
    ZmwType ScrapZmwType(void) const;

    /// \returns this record's average signal-to-noise for each of A, C, G,
    ///          and T
    ///
    std::vector<float> SignalToNoise(void) const;

    /// \returns this record's type
    /// \sa RecordType
    RecordType Type(void) const;

    /// \}

public:
    /// \name Mapping Data
    /// \{

    /// \returns the record's aligned end position
    ///
    /// \note AlignedEnd is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Position AlignedEnd(void) const;

    /// \returns the record's aligned start position
    ///
    /// \note AlignedStart is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Position AlignedStart(void) const;

    /// \returns the record's strand as a Strand enum value
    Strand AlignedStrand(void) const;

    /// \returns the record's CIGAR data as a Cigar object
    ///
    /// \param[in] exciseAllClips   if true, remove all clipping operations
    ///                             (hard & soft) [default:false]
    ///
    Cigar CigarData(bool exciseAllClips = false) const;

    /// \returns true if this record was mapped by aligner
    bool IsMapped(void) const;

    /// \returns this record's mapping quality. A value of 255 indicates
    ///          "unknown"
    ///
    uint8_t MapQuality(void) const;

    /// \returns the number of deleted bases (relative to reference)
    size_t NumDeletedBases(void) const;

    /// \returns the number of inserted bases (relative to reference)
    size_t NumInsertedBases(void) const;

    /// \returns the number of matching bases (sum of '=' CIGAR op lengths)
    size_t NumMatches(void) const;

    /// \returns a tuple containing NumMatches (first) and NumMismatches
    ///         (second)
    ///
    std::pair<size_t, size_t> NumMatchesAndMismatches(void) const;

    /// \returns the number of mismatching bases (sum of 'X' CIGAR op lengths)
    size_t NumMismatches(void) const;

    /// \returns this record's reference ID, or -1 if unmapped.
    ///
    /// \note This is only a valid identifier within this %BAM file
    ///
    int32_t ReferenceId(void) const;

    /// \returns this record's reference name.
    ///
    /// \throws an exception if unmapped record.
    ///
    std::string ReferenceName(void) const;

    /// \returns the record's reference end position, or UnmappedPosition if
    ///          unmapped
    ///
    /// \note ReferenceEnd is in reference coordinates, NOT polymerase read
    ///       coordinates.
    ///
    Position ReferenceEnd(void) const;

    /// \returns the record's reference start position, or UnmappedPosition if
    ///          unmapped
    ///
    /// \note ReferenceStart is in reference coordinates, NOT polymerase read
    ///       coordinates.
    ///
    Position ReferenceStart(void) const;

    /// \}

public:
    /// \name Barcode Data
    /// \{

    /// \returns forward barcode id
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    int16_t BarcodeForward(void) const;

    /// \returns barcode call confidence (Phred-scaled posterior probability
    ///          of correct barcode call)
    ///
    /// \sa HasBarcodeQuality
    ///
    uint8_t BarcodeQuality(void) const;

    /// \returns reverse barcode id
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    int16_t BarcodeReverse(void) const;

    /// \returns the forward and reverse barcode ids
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    std::pair<int16_t,int16_t> Barcodes(void) const;

    /// \}

public:
    /// \name Auxiliary Data Queries
    /// \{

    /// \returns true if this record has AltLabelQV data
    bool HasAltLabelQV(void) const;

    /// \returns true if this record has AltLabelTag data
    bool HasAltLabelTag(void) const;

    /// \returns true if this record has Barcode data
    bool HasBarcodes(void) const;

    /// \returns true is this record has BarcodeQuality data
    bool HasBarcodeQuality(void) const;

    /// \returns true if this record has DeletionQV data
    bool HasDeletionQV(void) const;

    /// \returns true if this record has DeletionTag data
    bool HasDeletionTag(void) const;

    /// \returns true if this record has a HoleNumber
    bool HasHoleNumber(void) const;

    /// \returns true if this record has InsertionQV data
    bool HasInsertionQV(void) const;

    /// \returns true if this record has IPD data
    bool HasIPD(void) const;

    /// \returns true if this record has LabelQV data
    bool HasLabelQV(void) const;

    /// \returns true if this record has LocalContextFlags (absent in CCS)
    bool HasLocalContextFlags(void) const;

    /// \returns true if this record has MergeQV data
    bool HasMergeQV(void) const;

    /// \returns true if this record has NumPasses data
    bool HasNumPasses(void) const;

    /// \returns true if this record has Pkmean data
    bool HasPkmean(void) const;

    /// \returns true if this record has Pkmid data
    bool HasPkmid(void) const;

    /// \returns true if this record has Pkmean2 data
    bool HasPkmean2(void) const;

    /// \returns true if this record has Pkmid2 data
    bool HasPkmid2(void) const;

    /// \returns true if this record has PreBaseFrames aka IPD data
    bool HasPreBaseFrames(void) const;

    /// \returns true if this record has PrePulseFrames data
    bool HasPrePulseFrames(void) const;

    /// \returns true if this record has PulseCall data
    bool HasPulseCall(void) const;

    /// \returns true if this record has PulseCallWidth data
    bool HasPulseCallWidth(void) const;

    /// \returns true if this record has PulseMergeQV data
    bool HasPulseMergeQV(void) const;

    /// \returns true if this record has PulseWidth data
    bool HasPulseWidth(void) const;

    /// \returns true if this record has ReadAccuracyTag data
    bool HasReadAccuracy(void) const;

    /// \returns true if this record has QueryEnd data
    bool HasQueryEnd(void) const;

    /// \returns true if this record has QueryStart data
    bool HasQueryStart(void) const;

    /// \returns true if this record has ScrapRegionType data (only in SCRAP)
    bool HasScrapRegionType(void) const;

    /// \returns true if this record has scrap ZMW type data (only in SCRAP)
    bool HasScrapZmwType(void) const;

    /// \returns true if this record has signal-to-noise data (absent in
    ///          POLYMERASE)
    ///
    bool HasSignalToNoise(void) const;

    /// \returns true if this record has StartFrame data
    bool HasStartFrame(void) const;

    /// \returns true if this record has SubstitutionQV data
    bool HasSubstitutionQV(void) const;

    /// \returns true if this record has SubstitutionTag data
    bool HasSubstitutionTag(void) const;

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
    std::string AltLabelTag(Orientation orientation = Orientation::NATIVE) const;

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
    std::string DeletionTag(Orientation orientation = Orientation::NATIVE,
                            bool aligned = false,
                            bool exciseSoftClips = false) const;

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
    std::string Sequence(const Orientation orientation = Orientation::NATIVE,
                         bool aligned = false,
                         bool exciseSoftClips = false) const;

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
    std::string SubstitutionTag(Orientation orientation = Orientation::NATIVE,
                                bool aligned = false,
                                bool exciseSoftClips = false) const;

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
    QualityValues AltLabelQV(Orientation orientation = Orientation::NATIVE) const;

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
    QualityValues DeletionQV(Orientation orientation = Orientation::NATIVE,
                             bool aligned = false,
                             bool exciseSoftClips = false) const;

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
    QualityValues InsertionQV(Orientation orientation = Orientation::NATIVE,
                              bool aligned = false,
                              bool exciseSoftClips = false) const;

    /// \brief Fetches this record's LabelQV values ("pq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the
    ///       new QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    ///
    /// \returns LabelQV as QualityValues object
    ///
    QualityValues LabelQV(Orientation orientation = Orientation::NATIVE) const;

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
    QualityValues MergeQV(Orientation orientation = Orientation::NATIVE,
                          bool aligned = false,
                          bool exciseSoftClips = false) const;

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
    QualityValues Qualities(Orientation orientation = Orientation::NATIVE,
                            bool aligned = false,
                            bool exciseSoftClips = false) const;

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
    QualityValues SubstitutionQV(Orientation orientation = Orientation::NATIVE,
                                 bool aligned = false,
                                 bool exciseSoftClips = false) const;

    /// \}

public:
    /// \name Pulse Data
    /// \{

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
    Frames IPD(Orientation orientation = Orientation::NATIVE,
               bool aligned = false,
               bool exciseSoftClips = false) const;

    /// \brief Fetches this record's IPD values ("ip" tag), but does not upscale.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns IPD as Frames object
    ///
    Frames IPDRaw(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetches this record's Pkmean values ("pa" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmean as vector<float> object
    ///
    std::vector<float> Pkmean(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetches this record's Pkmid values ("pm" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmid as vector<float> object
    ///
    std::vector<float> Pkmid(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetches this record's Pkmean2 values ("pi" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmean as vector<float> object
    ///
    std::vector<float> Pkmean2(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetches this record's Pkmid2 values ("ps" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmid as vector<float> object
    ///
    std::vector<float> Pkmid2(Orientation orientation = Orientation::NATIVE) const;

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
    Frames PreBaseFrames(Orientation orientation = Orientation::NATIVE,
                         bool aligned = false,
                         bool exciseSoftClips = false) const;

    /// \brief Fetches this record's PrePulseFrames values ("pd" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PrePulseFrames as Frames object
    ///
    Frames PrePulseFrames(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetches this record's PulseCall values ("pc" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseCalls string
    ///
    std::string PulseCall(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetches this record's PulseCallWidth values ("px" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseCallWidth as Frames object
    ///
    Frames PulseCallWidth(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetch this record's PulseMergeQV values ("pg" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseMergeQV as QualityValues object
    ///
    QualityValues PulseMergeQV(Orientation orientation = Orientation::NATIVE) const;

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
    Frames PulseWidth(Orientation orientation = Orientation::NATIVE,
                      bool aligned = false,
                      bool exciseSoftClips = false) const;

    /// \brief Fetches this record's PulseWidth values ("pw" tag), but does not
    ///        upscale.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseWidth as Frames object
    ///
    Frames PulseWidthRaw(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Fetches this record's StartFrame values ("sf" tag).
    ///
    /// \param[in] orientation     Orientation of output
    ///
    /// \returns StartFrame as uint32_t vector
    ///
    std::vector<uint32_t> StartFrame(Orientation orientation = Orientation::NATIVE) const;

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
    const BamRecordImpl& Impl(void) const;

    /// \warning This method should be considered temporary and avoided as much
    ///          as possible. Direct access to the internal object is likely to
    ///          disappear as BamRecord interface matures.
    ///
    /// \returns reference to underlying BamRecordImpl object
    ///
    BamRecordImpl& Impl(void);

    /// \}

public:
    /// \name General Data
    /// \{

    /// \brief Sets this record's ZMW hole number.
    ///
    /// \param[in] holeNumber
    /// \returns reference to this record
    ///
    BamRecord& HoleNumber(const int32_t holeNumber);

    /// \brief Sets this record's local context flags
    ///
    /// \param[in] flags
    /// \returns reference to this record
    ///
    BamRecord& LocalContextFlags(const PacBio::BAM::LocalContextFlags flags);

    /// \brief Sets this record's "number of complete passes of the insert".
    ///
    /// \param[in] numPasses
    /// \returns reference to this record
    ///
    BamRecord& NumPasses(const int32_t numPasses);

    /// \brief Sets this record's query end position.
    ///
    /// \note Changing this will modify the name of non-CCS records.
    ///
    /// \param[in] pos
    /// \returns reference to this record
    ///
    BamRecord& QueryEnd(const PacBio::BAM::Position pos);

    /// \brief Sets this record's query start position.
    ///
    /// \note Changing this will modify the name of non-CCS records.
    ///
    /// \param[in] pos
    /// \returns reference to this record
    ///
    BamRecord& QueryStart(const PacBio::BAM::Position pos);

    /// \brief Sets this record's expected read accuracy [0, 1000]
    ///
    /// \param[in] accuracy
    /// \returns reference to this record
    ///
    BamRecord& ReadAccuracy(const Accuracy& accuracy);

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
    BamRecord& ScrapRegionType(const VirtualRegionType type);

    /// \brief Sets this scrap record's ScrapRegionType
    ///
    /// \param[in] type character equivalent of VirtualRegionType
    /// \returns reference to this record
    ///
    BamRecord& ScrapRegionType(const char type);

    /// \brief Sets this scrap record's ScrapZmwType
    ///
    /// \param[in] type
    /// \returns reference to this record
    ///
    BamRecord& ScrapZmwType(const ZmwType type);

    /// \brief Sets this scrap record's ScrapZmwType
    ///
    /// \param[in] type character equivalent of ZmwType
    /// \returns reference to this record
    ///
    BamRecord& ScrapZmwType(const char type);

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
    BamRecord& Barcodes(const std::pair<int16_t, int16_t>& barcodeIds);

    /// \brief Sets this record's barcode quality ('bq' tag)
    ///
    /// \param[in] quality Phred-scaled confidence call
    /// \returns reference to this record
    ///
    BamRecord& BarcodeQuality(const uint8_t quality);

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
    BamRecord& AltLabelQV(const QualityValues& altLabelQVs);

    /// \brief Sets this record's DeletionQV values ("dq" tag).
    ///
    /// \param[in] deletionQVs
    /// \returns reference to this record
    ///
    BamRecord& DeletionQV(const QualityValues& deletionQVs);

    /// \brief Sets this record's InsertionQV values ("iq" tag).
    ///
    /// \param[in] insertionQVs
    /// \returns reference to this record
    ///
    BamRecord& InsertionQV(const QualityValues& insertionQVs);

    /// \brief Sets this record's LabelQV values ("pq" tag).
    ///
    /// \param[in] labelQVs
    /// \returns reference to this record
    ///
    BamRecord& LabelQV(const QualityValues& labelQVs);

    /// \brief Sets this record's MergeQV values ("mq" tag).
    ///
    /// \param[in] mergeQVs
    /// \returns reference to this record
    ///
    BamRecord& MergeQV(const QualityValues& mergeQVs);

    /// \brief Sets this record's SubstitutionQV values ("sq" tag).
    ///
    /// \param[in] substitutionQVs
    /// \returns reference to this record
    ///
    BamRecord& SubstitutionQV(const QualityValues& substitutionQVs);

    /// \}

public:
    /// \name Pulse Data
    /// \{

    /// \brief Sets this record's IPD values ("ip" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& IPD(const Frames& frames,
                   const FrameEncodingType encoding);

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
    BamRecord& Pkmean(const std::vector<uint16_t>& encodedPhotons);

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
    BamRecord& Pkmid(const std::vector<uint16_t>& encodedPhotons);

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
    BamRecord& Pkmean2(const std::vector<uint16_t>& encodedPhotons);

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
    BamRecord& Pkmid2(const std::vector<uint16_t>& encodedPhotons);

    /// \brief Sets this record's PreBaseFrames aka IPD values ("ip" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PreBaseFrames(const Frames& frames,
                             const FrameEncodingType encoding);

    /// \brief Sets this record's PrePulseFrames values ("pd" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PrePulseFrames(const Frames& frames,
                              const FrameEncodingType encoding);

    /// \brief Sets this record's PulseCall values ("pc" tag).
    ///
    /// \param[in] tags
    /// \returns reference to this record
    ///
    BamRecord& PulseCall(const std::string& tags);

    /// \brief Sets this record's PulseCallWidth values ("px" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PulseCallWidth(const Frames& frames,
                              const FrameEncodingType encoding);

    /// \brief Sets this record's PulseMergeQV values ("pg" tag).
    ///
    /// \param[in] pulseMergeQVs
    /// \returns reference to this record
    ///
    BamRecord& PulseMergeQV(const QualityValues& pulseMergeQVs);

    /// \brief Sets this record's PulseWidth values ("pw" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or
    ///                     16-bit lossless)
    /// \returns reference to this record
    ///
    BamRecord& PulseWidth(const Frames& frames,
                          const FrameEncodingType encoding);

    /// \brief Sets this record's StartFrame values ("sf" tag).
    ///
    /// \param[in] startFrame
    /// \returns reference to this record
    ///
    BamRecord& StartFrame(const std::vector<uint32_t>& startFrame);

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
    void ResetCachedPositions(void) const;

    /// \brief Resets cached aligned start/end.
    ///
    /// \note This method should not be needed in most client code. It exists
    ///       primarily as a hook for internal reading loops (queries, index
    ///       build, etc.) It's essentially a workaround and will likely be
    ///       removed from the API.
    ///
    void ResetCachedPositions(void);

    /// \brief Updates the record's name (BamRecord::FullName) to reflect
    ///        modifications to name components (movie name, ZMW hole number,
    ///        etc.)
    ///
    void UpdateName(void);

    /// \}

public:
    /// \name Pulse Data
    /// \{

    static const float photonFactor;

    static std::vector<uint16_t> EncodePhotons(const std::vector<float>& data);

    /// \}

public:
    /// \name Clipping & Mapping
    /// \{

    /// Creates a copied record from input, with clipping applied
    static BamRecord Clipped(const BamRecord& input,
                             const ClipType clipType,
                             const PacBio::BAM::Position start,
                             const PacBio::BAM::Position end);

    /// Creates a copied record from input, with mapping applied
    static BamRecord Mapped(const BamRecord& input,
                            const int32_t referenceId,
                            const Position refStart,
                            const Strand strand,
                            const Cigar& cigar,
                            const uint8_t mappingQuality);

    /// Applies clipping to this record
    BamRecord& Clip(const ClipType clipType,
                    const PacBio::BAM::Position start,
                    const PacBio::BAM::Position end);

    /// Creates a copied record from this one, with clipping applied
    BamRecord Clipped(const ClipType clipType,
                      const PacBio::BAM::Position start,
                      const PacBio::BAM::Position end) const;

    /// Applies mapping to this record
    BamRecord& Map(const int32_t referenceId,
                   const Position refStart,
                   const Strand strand,
                   const Cigar& cigar,
                   const uint8_t mappingQuality);

    /// Creates a copied record from this one, with mapping applied
    BamRecord Mapped(const int32_t referenceId,
                     const Position refStart,
                     const Strand strand,
                     const Cigar& cigar,
                     const uint8_t mappingQuality) const;
    /// \}

private:
    BamRecordImpl impl_;

public:
    /// public & mutable so that queries can directly set the header info,
    /// even on a record that is const from client code's perspective
    mutable BamHeader header_;

private:
    /// \internal
    /// cached positions (mutable to allow lazy-calc in const methods)
    mutable Position alignedStart_;
    mutable Position alignedEnd_;

private:
    /// \internal
    std::vector<float> FetchPhotons(const std::string& tagName,
                                    const Orientation orientation) const;
    std::string FetchBasesRaw(const std::string& tagName) const;

    std::string FetchBases(const std::string& tagName,
                           const Orientation orientation) const;

    std::string FetchBases(const std::string& tagName,
                           const Orientation orientation,
                           const bool aligned,
                           const bool exciseSoftClips) const;

    Frames FetchFramesRaw(const std::string& tagName) const;

    Frames FetchFrames(const std::string& tagName,
                       const Orientation orientation) const;

    Frames FetchFrames(const std::string& tagName,
                       const Orientation orientation,
                       const bool aligned,
                       const bool exciseSoftClips) const;

    QualityValues FetchQualitiesRaw(const std::string& tagName) const;

    QualityValues FetchQualities(const std::string& tagName,
                                 const Orientation orientation) const;

    QualityValues FetchQualities(const std::string& tagName,
                                 const Orientation orientation,
                                 const bool aligned,
                                 const bool exciseSoftClips) const;

    void ClipFields(const size_t clipPos, const size_t clipLength);
    BamRecord& ClipToQuery(const PacBio::BAM::Position start,
                           const PacBio::BAM::Position end);
    BamRecord& ClipToReference(const PacBio::BAM::Position start,
                               const PacBio::BAM::Position end);
    BamRecord& ClipToReferenceForward(const PacBio::BAM::Position start,
                                      const PacBio::BAM::Position end);
    BamRecord& ClipToReferenceReverse(const PacBio::BAM::Position start,
                                      const PacBio::BAM::Position end);

private:
    // marked const to allow calling from const methods
    // but updates our mutable cached values
    void CalculateAlignedPositions(void) const;

    friend class internal::BamRecordMemory;
};

/// \brief Provides a re-usable "view" onto a BamRecord
///
/// This class acts a convenience wrapper for working with per-base BamRecord
/// data. Most of these BamRecord methods take a list of parameters, to adjust
/// how the underlying data are presented to client code. Often these parameters
/// will be re-used for each BamRecord method call. Thus, to simplify such
/// client code, a BamRecordView can be used to state those parameters once, and
/// then simply request the desired fields.
///
/// \internal
/// \todo Sync up method names with BamRecord
/// \endinternal
///
class PBBAM_EXPORT BamRecordView
{
public:
    /// \brief Constructs a view onto \p record using the supplied parameters.
    ///
    /// For frame or QV data, if \p aligned is true, a value of 0 (Accuracy or
    /// QualityValue) will be used at each inserted or padded base location.
    ///
    /// \param[in] record           BamRecord data source.
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    BamRecordView(const BamRecord& record,
                  const Orientation orientation,
                  const bool aligned,
                  const bool exciseSoftClips);

public:

    /// \returns BamRecord::AltLabelQV with this view's parameters applied
    QualityValues AltLabelQVs(void) const;

    /// \returns BamRecord::AltLabelTag with this view's parameters applied
    std::string AltLabelTags(void) const;

    /// \returns BamRecord::DeletionQV with this view's parameters applied
    QualityValues DeletionQVs(void) const;

    /// \returns BamRecord::DeletionTag with this view's parameters applied
    std::string DeletionTags(void) const;

    /// \returns BamRecord::InsertionQV with this view's parameters applied
    QualityValues InsertionQVs(void) const;

    /// \returns BamRecord::IPD with this view's parameters applied
    Frames IPD(void) const;

    /// \returns BamRecord::LabelQV with this view's parameters applied
    QualityValues LabelQVs(void) const;

    /// \returns BamRecord::MergeQV with this view's parameters applied
    QualityValues MergeQVs(void) const;

    /// \returns BamRecord::PulseMergeQV with this view's parameters applied
    QualityValues PulseMergeQVs(void) const;

    /// \returns BamRecord::Pkmean with this view's parameters applied
    std::vector<float> Pkmean(void) const;

    /// \returns BamRecord::Pkmid with this view's parameters applied
    std::vector<float> Pkmid(void) const;

    /// \returns BamRecord::Pkmean2 with this view's parameters applied
    std::vector<float> Pkmean2(void) const;

    /// \returns BamRecord::Pkmid2 with this view's parameters applied
    std::vector<float> Pkmid2(void) const;

    /// \returns BamRecord::PreBaseFrames with this view's parameters applied
    Frames PrebaseFrames(void) const;

    /// \returns BamRecord::PrePulseFrames with this view's parameters applied
    Frames PrePulseFrames(void) const;

    /// \returns BamRecord::PulseCalls with this view's parameters applied
    std::string PulseCalls(void) const;

    /// \returns BamRecord::PulseCallWidth with this view's parameters applied
    Frames PulseCallWidth(void) const;

    /// \returns BamRecord::PulseWidths with this view's parameters applied
    Frames PulseWidths(void) const;

    /// \returns BamRecord::Qualities with this view's parameters applied
    QualityValues Qualities(void) const;

    /// \returns BamRecord::Sequence with this view's parameters applied
    std::string Sequence(void) const;

    /// \returns BamRecord::StartFrame with this view's parameters applied
    std::vector<uint32_t> StartFrames(void) const;

    /// \returns BamRecord::SubstitutionQV with this view's parameters applied
    QualityValues SubstitutionQVs(void) const;

    /// \returns BamRecord::SubstitutionTag with this view's parameters applied
    std::string SubstitutionTags(void) const;

private:
    const BamRecord& record_;
    Orientation orientation_;
    bool aligned_;
    bool exciseSoftClips_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/BamRecord.inl"

#endif // BAMRECORD_H
