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
/// \file BamRecord.inl
/// \brief Inline implementations for the BamRecord & BamRecordView classes.
//
// Author: Derek Barnett

#include "pbbam/BamRecord.h"

namespace PacBio {
namespace BAM {

inline BamRecord BamRecord::Clipped(const BamRecord& input,
                                    const ClipType clipType,
                                    const PacBio::BAM::Position start,
                                    const PacBio::BAM::Position end)
{
    return input.Clipped(clipType, start, end);
}

inline BamRecord BamRecord::Clipped(const ClipType clipType,
                                    const PacBio::BAM::Position start,
                                    const PacBio::BAM::Position end) const
{
    BamRecord result(*this);
    result.Clip(clipType, start, end);
    return result;
}

inline BamRecord BamRecord::Mapped(const BamRecord& input,
                                   const int32_t referenceId,
                                   const Position refStart,
                                   const Strand strand,
                                   const Cigar& cigar,
                                   const uint8_t mappingQuality)
{
    return input.Mapped(referenceId, refStart, strand, cigar, mappingQuality);
}

inline BamRecord BamRecord::Mapped(const int32_t referenceId,
                                   const Position refStart,
                                   const Strand strand,
                                   const Cigar& cigar,
                                   const uint8_t mappingQuality) const
{
    BamRecord result(*this);
    result.Map(referenceId, refStart, strand, cigar, mappingQuality);
    return result;
}


inline BamRecordView::BamRecordView(const BamRecord& record,
                                    const Orientation orientation,
                                    const bool aligned,
                                    const bool exciseSoftClips)
    : record_(record)
    , orientation_(orientation)
    , aligned_(aligned)
    , exciseSoftClips_(exciseSoftClips)
{ }

inline QualityValues BamRecordView::AltLabelQVs(void) const
{ return record_.AltLabelQV(orientation_); }

inline std::string BamRecordView::AltLabelTags(void) const
{ return record_.AltLabelTag(orientation_); }

inline QualityValues BamRecordView::DeletionQVs(void) const
{ return record_.DeletionQV(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::DeletionTags(void) const
{ return record_.DeletionTag(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::InsertionQVs(void) const
{ return record_.InsertionQV(orientation_, aligned_, exciseSoftClips_); }

inline Frames BamRecordView::IPD(void) const
{ return record_.IPD(orientation_, aligned_, exciseSoftClips_); }

inline Frames BamRecordView::PrebaseFrames(void) const
{ return record_.IPD(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::LabelQVs(void) const
{ return record_.LabelQV(orientation_); }

inline QualityValues BamRecordView::MergeQVs(void) const
{ return record_.MergeQV(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::PulseMergeQVs(void) const
{ return record_.PulseMergeQV(orientation_); }

inline std::vector<float> BamRecordView::Pkmean(void) const
{ return record_.Pkmean(orientation_); }

inline std::vector<float> BamRecordView::Pkmid(void) const
{ return record_.Pkmid(orientation_); }

inline Frames BamRecordView::PrePulseFrames(void) const
{ return record_.PrePulseFrames(orientation_); }

inline std::string BamRecordView::PulseCalls(void) const
{ return record_.PulseCall(orientation_); }

inline Frames BamRecordView::PulseCallWidth(void) const
{ return record_.PulseCallWidth(orientation_); }

inline Frames BamRecordView::PulseWidths(void) const
{ return record_.PulseWidth(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::Qualities(void) const
{ return record_.Qualities(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::Sequence(void) const
{ return record_.Sequence(orientation_, aligned_, exciseSoftClips_); }

inline std::vector<uint32_t> BamRecordView::StartFrames(void) const
{ return record_.StartFrame(orientation_); }

inline QualityValues BamRecordView::SubstitutionQVs(void) const
{ return record_.SubstitutionQV(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::SubstitutionTags(void) const
{ return record_.SubstitutionTag(orientation_, aligned_, exciseSoftClips_); }

} // namespace BAM
} // namespace PacBio
