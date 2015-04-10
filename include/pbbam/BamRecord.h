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

// Author: Derek Barnett

#ifndef BAMRECORD_H
#define BAMRECORD_H

#include "pbbam/Accuracy.h"
#include "pbbam/Frames.h"
#include "pbbam/BamRecordImpl.h"
#include "pbbam/BamHeader.h"
#include "pbbam/Orientation.h"
#include "pbbam/ReadGroupInfo.h"
#include "pbbam/Strand.h"
#include "pbbam/QualityValues.h"
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal { class BamRecordMemory; }

enum class ClipType
{
    CLIP_NONE
  , CLIP_TO_QUERY
  , CLIP_TO_REFERENCE
};

enum class RecordType
{
    POLYMERASE
  , HQREGION
  , SUBREAD
  , CCS
  , SCRAP
  , UNKNOWN
};

enum class FrameEncodingType
{
    LOSSY
  , LOSSLESS
};

class PBBAM_EXPORT BamRecord
{
public:
    /// \name Constructors & Related Methods
    /// \{

    BamRecord(void);
    BamRecord(const BamHeader::SharedPtr& header);
    BamRecord(const BamRecordImpl& impl);
    BamRecord(BamRecordImpl&& impl);
    BamRecord(const BamRecord& other);
    BamRecord(BamRecord&& other);
    BamRecord& operator=(const BamRecord& other);
    BamRecord& operator=(BamRecord&& other);
    virtual ~BamRecord(void);

    /// \}

public:
    /// \name Per-Record Data
    /// \{

    /// \note AlignedStart is in polymerase read coordinates, NOT genomic coordinates.
    ///
    /// \returns the record's aligned start position
    Position AlignedStart(void) const;

    /// \note AlignedEnd is in polymerase read coordinates, NOT genomic coordinates.
    ///
    /// \returns the record's aligned end position
    Position AlignedEnd(void) const;

    /// \returns the record's strand as a Strand enum value
    Strand AlignedStrand(void) const;

    /// \returns the record's CIGAR data as a Cigar object
    Cigar CigarData(void) const;

    /// \returns this record's full name
    /// \sa BamRecordImpl::Name
    std::string FullName(void) const;

    /// \returns true if this record has DeletionQV data
    bool HasDeletionQV(void) const;

    /// \returns true if this record has DeletionTag data
    bool HasDeletionTag(void) const;

    /// \returns true if this record has InsertionQV data
    bool HasInsertionQV(void) const;

    /// \returns true if this record has IPD data
    bool HasIPD(void) const;

    /// \returns true if this record has MergeQV data
    bool HasMergeQV(void) const;

    /// \returns true if this record has PulseWidth data
    bool HasPulseWidth(void) const;

    /// \returns true if this record has SubstitutionQV data
    bool HasSubstitutionQV(void) const;

    /// \returns true if this record has SubstitutionTag data
    bool HasSubstitutionTag(void) const;

    /// \returns shared pointer to this record's associated BamHeader
    BamHeader::SharedPtr Header(void) const;

    /// \returns ZMW hole number
    int32_t HoleNumber(void) const;

    /// \returns true if this record was mapped by aligner
    /// \sa BamRecordImpl::IsMapped
    bool IsMapped(void) const;

    /// \returns this record's mapping quality. A value of 255 indicates "unknown"
    uint8_t MapQuality(void) const;

    /// \returns this record's movie name
    std::string MovieName(void) const;

    /// \returns "number of complete passes of the insert"
    int32_t NumPasses(void) const;

    /// \note QueryStart is in polymerase read coordinates, NOT genomic coordinates.
    ///
    /// \returns the record's query start position
    Position QueryStart(void) const;

    /// \note QueryEnd is in polymerase read coordinates, NOT genomic coordinates.
    ///
    /// \returns the record's query end position
    Position QueryEnd(void) const;

    /// \returns this record's expected read accuracy [0, 1000]
    Accuracy ReadAccuracy(void) const;

    /// \returns ReadGroupInfo object for this record
    ReadGroupInfo ReadGroup(void) const;

    /// \returns ID of this record's read group
    /// \sa ReadGroupInfo::Id
    std::string ReadGroupId(void) const;

    /// \returns this record's reference ID, or -1 if unmapped
    int32_t ReferenceId(void) const;

    /// \note ReferenceStart is in reference coordinates, NOT polymerase read coordinates.
    ///
    /// \returns the record's reference start position, or UnmappedPosition if unmapped
    Position ReferenceStart(void) const;

    /// \note ReferenceEnd is in reference coordinates, NOT polymerase read coordinates.
    ///
    /// \returns the record's reference end position, or UnmappedPosition if unmapped
    Position ReferenceEnd(void) const;

    /// \returns this record's type
    /// \sa RecordType
    RecordType Type(void) const;

    /// \}

public:
    /// \name Per-Base Data
    /// \{

    /// \brief Fetch this record's DeletionQV values ("dq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns DeletionQV as QualityValues object
    ///
    QualityValues DeletionQV(Orientation orientation = Orientation::NATIVE,
                             bool aligned = false,
                             bool exciseSoftClips = false) const;

    /// \brief Fetch this record's DeletionTag values ("dt" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       gap chars will be '-' and padding chars will be '*'.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns DeletionTag string
    ///
    std::string DeletionTag(Orientation orientation = Orientation::NATIVE,
                            bool aligned = false,
                            bool exciseSoftClips = false) const;

    /// \brief Fetch this record's InsertionQV values ("iq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns InsertionQVs as QualityValues object
    ///
    QualityValues InsertionQV(Orientation orientation = Orientation::NATIVE,
                              bool aligned = false,
                              bool exciseSoftClips = false) const;

    /// \brief Fetch this record's IPD values ("ip" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       frames will have a value of 0;
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns IPD as Frames object
    ///
    Frames IPD(Orientation orientation = Orientation::NATIVE,
               bool aligned = false,
               bool exciseSoftClips = false) const;

    /// \brief Fetch this record's MergeQV values ("mq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns MergeQV as QualityValues object
    ///
    QualityValues MergeQV(Orientation orientation = Orientation::NATIVE,
                          bool aligned = false,
                          bool exciseSoftClips = false) const;

    /// \brief Fetch this record's PulseWidth values ("pw" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       frames will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns PulseWidths as Frames object
    ///
    Frames PulseWidth(Orientation orientation = Orientation::NATIVE,
                      bool aligned = false,
                      bool exciseSoftClips = false) const;

    /// \brief Fetch this record's BAM quality values (QUAL field).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns BAM qualities as QualityValues object
    ///
    QualityValues Qualities(Orientation orientation = Orientation::NATIVE,
                            bool aligned = false,
                            bool exciseSoftClips = false) const;

    /// \brief Fetch this record's DNA sequence (SEQ field).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       gap chars will be '-' and padding chars will be '*'.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns sequence string
    ///
    std::string Sequence(const Orientation orientation = Orientation::NATIVE,
                         bool aligned = false,
                         bool exciseSoftClips = false) const;

    /// \brief Fetch this record's SubstitutionQV values ("sq" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       QVs will have a value of 0.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns SubstitutionQV as QualityValues object
    ///
    QualityValues SubstitutionQV(Orientation orientation = Orientation::NATIVE,
                                 bool aligned = false,
                                 bool exciseSoftClips = false) const;

    /// \brief Fetch this record's SubstitutionTag values ("st" tag).
    ///
    /// \note If \p aligned is true, and gaps/padding need to be inserted, the new
    ///       gap chars will be '-' and padding chars will be '*'.
    ///
    /// \param[in] orientation     Orientation of output.
    /// \param[in] aligned         if true, gaps/padding will be inserted, per Cigar info.
    /// \param[in] exciseSoftClips if true, any soft-clipped positions will be removed from query ends
    ///
    /// \returns SubstitutionTags string
    ///
    std::string SubstitutionTag(Orientation orientation = Orientation::NATIVE,
                                bool aligned = false,
                                bool exciseSoftClips = false) const;

    /// \}

public:
    /// \name Low-Level
    /// \{

    /// \warning This method should be considered temporary and avoided as much as possible.
    ///          Direct access to the internal object is likely to disappear as BamRecord interface matures.
    ///
    /// \returns const reference to underlying BamRecordImpl object
    const BamRecordImpl& Impl(void) const;

    /// \warning This method should be considered temporary and avoided as much as possible.
    ///          Direct access to the internal object is likely to disappear as BamRecord interface matures.
    ///
    /// \returns reference to underlying BamRecordImpl object
    BamRecordImpl& Impl(void);

    /// \}

public:
    /// \name Per-Record Data
    /// \{
    ///

    /// Sets this record's ZMW hole number.
    ///
    /// \param[in] numPasses
    /// \returns reference to this record
    BamRecord& HoleNumber(const int32_t holeNumber);

    /// Sets this record's "number of complete passes of the insert".
    ///
    /// \param[in] numPasses
    /// \returns reference to this record
    BamRecord& NumPasses(const int32_t numPasses);

    /// Sets this record's expected read accuracy [0, 1000]
    ///
    /// \param[in] accuracy
    /// \returns reference to this record
    BamRecord& ReadAccuracy(const Accuracy& accuracy);

    /// \}

public:
    /// \name Per-Base Data
    /// \{

    /// Sets this record's DeletionQV values ("dq" tag).
    ///
    /// \param[in] deletionQVs
    /// \returns reference to this record
    BamRecord& DeletionQV(const QualityValues& deletionQVs);

    /// Sets this record's DeletionTag values ("dt" tag).
    ///
    /// \param[in] tags
    /// \returns reference to this record
    BamRecord& DeletionTag(const std::string& tags);

    /// Sets this record's InsertionQV values ("iq" tag).
    ///
    /// \param[in] insertionQVs
    /// \returns reference to this record
    BamRecord& InsertionQV(const QualityValues& insertionQVs);

    /// Sets this record's IPD values ("ip" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or 16-bit lossless)
    /// \returns reference to this record
    BamRecord& IPD(const Frames& frames,
                   const FrameEncodingType encoding);

    /// Sets this record's MergeQV values ("mq" tag).
    ///
    /// \param[in] mergeQVs
    /// \returns reference to this record
    BamRecord& MergeQV(const QualityValues& mergeQVs);

    /// Sets this record's PulseWidth values ("pw" tag).
    ///
    /// \param[in] frames
    /// \param[in] encoding specify how to encode the data (8-bit lossy, or 16-bit lossless)
    /// \returns reference to this record
    BamRecord& PulseWidth(const Frames& frames,
                          const FrameEncodingType encoding);

    /// Sets this record's SubstitutionQV values ("sq" tag).
    ///
    /// \param[in] substitutionQVs
    /// \returns reference to this record
    BamRecord& SubstitutionQV(const QualityValues& substitutionQVs);

    /// Sets this record's SubstitutionTag values ("st" tag).
    ///
    /// \param[in] tags
    /// \returns reference to this record
    BamRecord& SubstitutionTag(const std::string& tags);

    /// \}


//public:
//    BamRecord& QueryEnd(const PacBio::BAM::Position pos);
//    BamRecord& QueryStart(const PacBio::BAM::Position pos);


//    BamRecord& MovieName(const std::string& movie);

//    BamRecord& ReadGroup(const ReadGroupInfo& rg);
//    BamRecord& ReadGroupId(const std::string& id);
//    BamRecord& ReferenceStart(const PacBio::BAM::Position pos);

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
    BamHeader::SharedPtr header_;

    // cached positions (mutable to allow lazy-calc in const methods)
    mutable Position alignedStart_;
    mutable Position alignedEnd_;

private:
    std::string FetchBases(const std::string& tagName,
                           const Orientation orientation,
                           const bool aligned,
                           const bool exciseSoftClips) const;

    Frames FetchFrames(const std::string& tagName,
                       const Orientation orientation,
                       const bool aligned,
                       const bool exciseSoftClips) const;

    QualityValues FetchQualities(const std::string& tagName,
                                 const Orientation orientation,
                                 const bool aligned,
                                 const bool exciseSoftClips) const;

private:
    // marked const to allow calling from const methods
    //   but lazy-calc mutable, cached values
    void CalculateAlignedPositions(void) const;

//    void UpdateName(void);

    friend class internal::BamRecordMemory;
};

inline
BamRecord BamRecord::Clipped(const BamRecord& input,
                             const ClipType clipType,
                             const PacBio::BAM::Position start,
                             const PacBio::BAM::Position end)
{
    return input.Clipped(clipType, start, end);
}

inline
BamRecord BamRecord::Clipped(const ClipType clipType,
                             const PacBio::BAM::Position start,
                             const PacBio::BAM::Position end) const
{
    BamRecord result(*this);
    result.Clip(clipType, start, end);
    return result;
}

inline
BamRecord BamRecord::Mapped(const BamRecord& input,
                            const int32_t referenceId,
                            const Position refStart,
                            const Strand strand,
                            const Cigar& cigar,
                            const uint8_t mappingQuality)
{
    return input.Mapped(referenceId, refStart, strand, cigar, mappingQuality);
}

inline
BamRecord BamRecord::Mapped(const int32_t referenceId,
                            const Position refStart,
                            const Strand strand,
                            const Cigar& cigar,
                            const uint8_t mappingQuality) const
{
    BamRecord result(*this);
    result.Map(referenceId, refStart, strand, cigar, mappingQuality);
    return result;
}

class PBBAM_EXPORT BamRecordView
{
public:
    BamRecordView(const BamRecord& record,
                  const Orientation orientation,
                  const bool aligned,
                  const bool exciseSoftClips)
        : record_(record)
        , orientation_(orientation)
        , aligned_(aligned)
        , exciseSoftClips_(exciseSoftClips)
    { }

public:
    QualityValues DeletionQVs(void) const
    { return record_.DeletionQV(orientation_, aligned_, exciseSoftClips_); }

    std::string DeletionTags(void) const
    { return record_.DeletionTag(orientation_, aligned_, exciseSoftClips_); }

    QualityValues InsertionQVs(void) const
    { return record_.InsertionQV(orientation_, aligned_, exciseSoftClips_); }

    Frames IPD(void) const
    { return record_.IPD(orientation_, aligned_, exciseSoftClips_); }

    QualityValues MergeQVs(void) const
    { return record_.MergeQV(orientation_, aligned_, exciseSoftClips_); }

    Frames PulseWidths(void) const
    { return record_.PulseWidth(orientation_, aligned_, exciseSoftClips_); }

    QualityValues Qualities(void) const
    { return record_.Qualities(orientation_, aligned_, exciseSoftClips_); }

    std::string Sequence(void) const
    { return record_.Sequence(orientation_, aligned_, exciseSoftClips_); }

    QualityValues SubstitutionQVs(void) const
    { return record_.SubstitutionQV(orientation_, aligned_, exciseSoftClips_); }

    std::string SubstitutionTags(void) const
    { return record_.SubstitutionTag(orientation_, aligned_, exciseSoftClips_); }

private:
    const BamRecord& record_;
    Orientation orientation_;
    bool aligned_;
    bool exciseSoftClips_;
};

} // namespace BAM
} // namespace PacBio

#endif // BAMRECORD_H
