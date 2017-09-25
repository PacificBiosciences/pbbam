// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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
/// \file BamRecordView.inl
/// \brief Inline implementations for the BamRecordView class.
//
// Author: Derek Barnett

#include "pbbam/BamRecordView.h"

namespace PacBio {
namespace BAM {

inline BamRecordView::BamRecordView(const BamRecord& record,
                                    const Orientation orientation,
                                    const bool aligned,
                                    const bool exciseSoftClips,
                                    const PulseBehavior pulseBehavior)
    : record_(record)
    , orientation_(orientation)
    , aligned_(aligned)
    , exciseSoftClips_(exciseSoftClips)
    , pulseBehavior_(pulseBehavior)
{ }

inline QualityValues BamRecordView::AltLabelQVs() const
{ return record_.AltLabelQV(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::string BamRecordView::AltLabelTags() const
{ return record_.AltLabelTag(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline QualityValues BamRecordView::DeletionQVs() const
{ return record_.DeletionQV(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::DeletionTags() const
{ return record_.DeletionTag(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::InsertionQVs() const
{ return record_.InsertionQV(orientation_, aligned_, exciseSoftClips_); }

inline Frames BamRecordView::IPD() const
{ return record_.IPD(orientation_, aligned_, exciseSoftClips_); }

inline Frames BamRecordView::PrebaseFrames() const
{ return record_.IPD(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::LabelQVs() const
{ return record_.LabelQV(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline QualityValues BamRecordView::MergeQVs() const
{ return record_.MergeQV(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::PulseMergeQVs() const
{ return record_.PulseMergeQV(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmean() const
{ return record_.Pkmean(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmid() const
{ return record_.Pkmid(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmean2() const
{ return record_.Pkmean2(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmid2() const
{ return record_.Pkmid2(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline Frames BamRecordView::PrePulseFrames() const
{ return record_.PrePulseFrames(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::string BamRecordView::PulseCalls() const
{ return record_.PulseCall(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline Frames BamRecordView::PulseCallWidth() const
{ return record_.PulseCallWidth(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline Frames BamRecordView::PulseWidths() const
{ return record_.PulseWidth(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::Qualities() const
{ return record_.Qualities(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::Sequence() const
{ return record_.Sequence(orientation_, aligned_, exciseSoftClips_); }

inline std::vector<uint32_t> BamRecordView::StartFrames() const
{ return record_.StartFrame(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline QualityValues BamRecordView::SubstitutionQVs() const
{ return record_.SubstitutionQV(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::SubstitutionTags() const
{ return record_.SubstitutionTag(orientation_, aligned_, exciseSoftClips_); }

} // namespace BAM
} // namespace PacBio
