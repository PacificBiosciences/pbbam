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
/// \brief Defines the BamRecord class.
//
// Author: Derek Barnett

#ifndef BAMRECORD_H
#define BAMRECORD_H

#include "pbbam/Accuracy.h"
#include "pbbam/Frames.h"
#include "pbbam/BamRecordImpl.h"
#include "pbbam/BamHeader.h"
#include "pbbam/ClipType.h"
#include "pbbam/FrameEncodingType.h"
#include "pbbam/LocalContextFlags.h"
#include "pbbam/Orientation.h"
#include "pbbam/PulseBehavior.h"
#include "pbbam/PulseExclusionReason.h"
#include "pbbam/ReadGroupInfo.h"
#include "pbbam/RecordType.h"
#include "pbbam/Strand.h"
#include "pbbam/QualityValues.h"
#include "pbbam/virtual/VirtualRegionType.h"
#include "pbbam/ZmwType.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal {

class BamRecordMemory;
class Pulse2BaseCache;

} // namespace internal

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
    BamRecord(BamRecord&& other);
    BamRecord& operator=(const BamRecord& other);
    BamRecord& operator=(BamRecord&& other);
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
    int32_t HoleNumber() const;

    /// \returns this record's LocalContextFlags
    PacBio::BAM::LocalContextFlags LocalContextFlags() const;

    /// \returns this record's movie name
    std::string MovieName() const;

    /// \returns "number of complete passes of the insert"
    int32_t NumPasses() const;

    /// \returns the record's query end position, or Sequence().length() if not
    ///          stored
    /// \note QueryEnd is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Position QueryEnd() const;

    /// \returns the record's query start position, or 0 if not stored
    ///
    /// \note QueryStart is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Position QueryStart() const;

    /// \returns this record's expected read accuracy [0, 1000]
    Accuracy ReadAccuracy() const;

    /// \returns ReadGroupInfo object for this record
    ReadGroupInfo ReadGroup() const;

    /// \returns string ID of this record's read group
    /// \sa ReadGroupInfo::Id
    ///
    std::string ReadGroupId() const;

    /// \returns integer value for this record's read group ID
    int32_t ReadGroupNumericId() const;

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
    Position AlignedEnd() const;

    /// \returns the record's aligned start position
    ///
    /// \note AlignedStart is in polymerase read coordinates, NOT genomic
    ///       coordinates.
    ///
    Position AlignedStart() const;

    /// \returns the record's strand as a Strand enum value
    Strand AlignedStrand() const;

    /// \returns the record's CIGAR data as a Cigar object
    ///
    /// \param[in] exciseAllClips   if true, remove all clipping operations
    ///                             (hard & soft) [default:false]
    ///
    Cigar CigarData(bool exciseAllClips = false) const;

    /// \returns true if this record was mapped by aligner
    bool IsMapped() const;

    /// \returns this record's mapping quality. A value of 255 indicates
    ///          "unknown"
    ///
    uint8_t MapQuality() const;

    /// \returns the number of deleted bases (relative to reference)
    size_t NumDeletedBases() const;

    /// \returns the number of inserted bases (relative to reference)
    size_t NumInsertedBases() const;

    /// \returns the number of matching bases (sum of '=' CIGAR op lengths)
    size_t NumMatches() const;

    /// \returns a tuple containing NumMatches (first) and NumMismatches
    ///         (second)
    ///
    std::pair<size_t, size_t> NumMatchesAndMismatches() const;

    /// \returns the number of mismatching bases (sum of 'X' CIGAR op lengths)
    size_t NumMismatches() const;

    /// \returns this record's reference ID, or -1 if unmapped.
    ///
    /// \note This is only a valid identifier within this %BAM file
    ///
    int32_t ReferenceId() const;

    /// \returns this record's reference name.
    ///
    /// \throws an exception if unmapped record.
    ///
    std::string ReferenceName() const;

    /// \returns the record's reference end position, or UnmappedPosition if
    ///          unmapped
    ///
    /// \note ReferenceEnd is in reference coordinates, NOT polymerase read
    ///       coordinates.
    ///
    Position ReferenceEnd() const;

    /// \returns the record's reference start position, or UnmappedPosition if
    ///          unmapped
    ///
    /// \note ReferenceStart is in reference coordinates, NOT polymerase read
    ///       coordinates.
    ///
    Position ReferenceStart() const;

    /// \}

public:
    /// \name Barcode Data
    /// \{

    /// \returns forward barcode id
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    int16_t BarcodeForward() const;

    /// \returns barcode call confidence (Phred-scaled posterior probability
    ///          of correct barcode call)
    ///
    /// \sa HasBarcodeQuality
    ///
    uint8_t BarcodeQuality() const;

    /// \returns reverse barcode id
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    int16_t BarcodeReverse() const;

    /// \returns the forward and reverse barcode ids
    ///
    /// \throws std::runtime_error if barcode data is absent or malformed.
    /// \sa HasBarcodes
    ///
    std::pair<int16_t,int16_t> Barcodes() const;

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
    bool HasPulseExclusion(void) const;

    /// \returns true if this record has PulseMergeQV data
    bool HasPulseMergeQV() const;

    /// \returns true if this record has PulseWidth data
    bool HasPulseWidth() const;

    /// \returns true if this record has ReadAccuracyTag data
    bool HasReadAccuracy() const;

    /// \returns true if this record has QueryEnd data
    bool HasQueryEnd() const;

    /// \returns true if this record has QueryStart data
    bool HasQueryStart() const;

    /// \returns true if this record has ScrapRegionType data (only in SCRAP)
    bool HasScrapRegionType() const;

    /// \returns true if this record has scrap ZMW type data (only in SCRAP)
    bool HasScrapZmwType() const;

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
    std::string AltLabelTag(Orientation orientation = Orientation::NATIVE,
                            bool aligned = false,
                            bool exciseSoftClips = false,
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
    QualityValues AltLabelQV(Orientation orientation = Orientation::NATIVE,
                             bool aligned = false,
                             bool exciseSoftClips = false,
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
    QualityValues LabelQV(Orientation orientation = Orientation::NATIVE,
                          bool aligned = false,
                          bool exciseSoftClips = false,
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
    std::vector<float> Pkmean(Orientation orientation = Orientation::NATIVE,
                              bool aligned = false,
                              bool exciseSoftClips = false,
                              PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's Pkmid values ("pm" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmid as vector<float> object
    ///
    std::vector<float> Pkmid(Orientation orientation = Orientation::NATIVE,
                             bool aligned = false,
                             bool exciseSoftClips = false,
                             PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's Pkmean2 values ("pi" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmean as vector<float> object
    ///
    std::vector<float> Pkmean2(Orientation orientation = Orientation::NATIVE,
                               bool aligned = false,
                               bool exciseSoftClips = false,
                               PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's Pkmid2 values ("ps" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns Pkmid as vector<float> object
    ///
    std::vector<float> Pkmid2(Orientation orientation = Orientation::NATIVE,
                              bool aligned = false,
                              bool exciseSoftClips = false,
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
    Frames PreBaseFrames(Orientation orientation = Orientation::NATIVE,
                         bool aligned = false,
                         bool exciseSoftClips = false) const;

    /// \brief Fetches this record's PrePulseFrames values ("pd" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PrePulseFrames as Frames object
    ///
    Frames PrePulseFrames(Orientation orientation = Orientation::NATIVE,
                          bool aligned = false,
                          bool exciseSoftClips = false,
                          PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PulseCall values ("pc" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseCalls string
    ///
    std::string PulseCall(Orientation orientation = Orientation::NATIVE,
                          bool aligned = false,
                          bool exciseSoftClips = false,
                          PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PulseCallWidth values ("px" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseCallWidth as Frames object
    ///
    Frames PulseCallWidth(Orientation orientation = Orientation::NATIVE,
                          bool aligned = false,
                          bool exciseSoftClips = false,
                          PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetches this record's PulseExclusionReason values ("pe" tag).
    ///
    /// \returns vector of pulse exclusion reason value
    ///
    std::vector<PacBio::BAM::PulseExclusionReason>
    PulseExclusionReason(Orientation orientation = Orientation::NATIVE,
                         bool aligned = false,
                         bool exciseSoftClips = false,
                         PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    /// \brief Fetch this record's PulseMergeQV values ("pg" tag).
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseMergeQV as QualityValues object
    ///
    QualityValues PulseMergeQV(Orientation orientation = Orientation::NATIVE,
                               bool aligned = false,
                               bool exciseSoftClips = false,
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
    Frames PulseWidth(Orientation orientation = Orientation::NATIVE,
                      bool aligned = false,
                      bool exciseSoftClips = false) const;

    /// \brief Fetches this record's PulseWidth values ("pw" tag), but does not
    ///        upscale.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \returns PulseWidth as Frames object
    ///
    Frames PulseWidthRaw(Orientation orientation = Orientation::NATIVE,
                         bool aligned = false,
                         bool exciseSoftClips = false) const;

    /// \brief Fetches this record's StartFrame values ("sf" tag).
    ///
    /// \param[in] orientation     Orientation of output
    ///
    /// \returns StartFrame as uint32_t vector
    ///
    std::vector<uint32_t> StartFrame(Orientation orientation = Orientation::NATIVE,
                                     bool aligned = false,
                                     bool exciseSoftClips = false,
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

    ///
    /// \\brief Sets this record's PulseExclusionReason values ("pe" tag).
    /// \param[in] reasons
    /// \return reference to this record
    ///
    BamRecord& PulseExclusionReason(const std::vector<PacBio::BAM::PulseExclusionReason>& reasons);

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
    /// pulse to bam mapping cache
    mutable std::unique_ptr<internal::Pulse2BaseCache> p2bCache_;

public:
    /// clips the PacBio tags to a specified length
    void ClipTags(const size_t clipPos, const size_t clipLength);

private:
    ///\internal
    /// clipping methods

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
    ///\internal
    /// raw tag data fetching

    // sequence tags
    std::string FetchBasesRaw(const BamRecordTag tag) const;
    std::string FetchBases(const BamRecordTag tag,
                           const Orientation orientation = Orientation::NATIVE,
                           const bool aligned = false,
                           const bool exciseSoftClips = false,
                           const PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // frame tags
    Frames FetchFramesRaw(const BamRecordTag tag) const;
    Frames FetchFrames(const BamRecordTag tag,
                       const Orientation orientation = Orientation::NATIVE,
                       const bool aligned = false,
                       const bool exciseSoftClips = false,
                       const PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // pulse tags
    std::vector<float> FetchPhotonsRaw(const BamRecordTag tag) const;
    std::vector<float> FetchPhotons(const BamRecordTag tag,
                                    const Orientation orientation = Orientation::NATIVE,
                                    const bool aligned = false,
                                    const bool exciseSoftClips = false,
                                    const PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // QV tags
    QualityValues FetchQualitiesRaw(const BamRecordTag tag) const;
    QualityValues FetchQualities(const BamRecordTag tag,
                                 const Orientation orientation = Orientation::NATIVE,
                                 const bool aligned = false,
                                 const bool exciseSoftClips = false,
                                 const PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // UInt tags (e.g. start frame)
    //
    // TODO (DB): clean this up w.r.t FetchUInt8s
    //
    std::vector<uint32_t> FetchUInt32sRaw(const BamRecordTag tag) const;
    std::vector<uint32_t> FetchUInt32s(const BamRecordTag tag,
                                       const Orientation orientation = Orientation::NATIVE,
                                       const bool aligned = false,
                                       const bool exciseSoftClips = false,
                                       const PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

    // UInt tags (e.g. pulse exclusion)
    //
    // ODO (DB): clean this up w.r.t FetchUInt32s
    //
    std::vector<uint8_t> FetchUInt8sRaw(const BamRecordTag tag) const;
    std::vector<uint8_t> FetchUInt8s(const BamRecordTag tag,
                                     const Orientation orientation = Orientation::NATIVE,
                                     const bool aligned = false,
                                     const bool exciseSoftClips = false,
                                     const PulseBehavior pulseBehavior = PulseBehavior::ALL) const;

private:
    ///\internal
    /// marked const to allow calling from const methods
    /// but updates our mutable cached values
    void CalculateAlignedPositions() const;
    void CalculatePulse2BaseCache() const;

    friend class internal::BamRecordMemory;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/BamRecord.inl"

#endif // BAMRECORD_H
