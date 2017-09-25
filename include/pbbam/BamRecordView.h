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
/// \file BamRecordView.h
/// \brief Defines the BamRecordView class.
//
// Author: Derek Barnett

#ifndef BAMRECORDVIEW_H
#define BAMRECORDVIEW_H

#include <cstdint>

#include "pbbam/BamRecord.h"

namespace PacBio {
namespace BAM {

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
                  const bool exciseSoftClips,
                  const PulseBehavior pulseBehavior = PulseBehavior::ALL);

public:

    /// \returns BamRecord::AltLabelQV with this view's parameters applied
    QualityValues AltLabelQVs() const;

    /// \returns BamRecord::AltLabelTag with this view's parameters applied
    std::string AltLabelTags() const;

    /// \returns BamRecord::DeletionQV with this view's parameters applied
    QualityValues DeletionQVs() const;

    /// \returns BamRecord::DeletionTag with this view's parameters applied
    std::string DeletionTags() const;

    /// \returns BamRecord::InsertionQV with this view's parameters applied
    QualityValues InsertionQVs() const;

    /// \returns BamRecord::IPD with this view's parameters applied
    Frames IPD() const;

    /// \returns BamRecord::LabelQV with this view's parameters applied
    QualityValues LabelQVs() const;

    /// \returns BamRecord::MergeQV with this view's parameters applied
    QualityValues MergeQVs() const;

    /// \returns BamRecord::PulseMergeQV with this view's parameters applied
    QualityValues PulseMergeQVs() const;

    /// \returns BamRecord::Pkmean with this view's parameters applied
    std::vector<float> Pkmean() const;

    /// \returns BamRecord::Pkmid with this view's parameters applied
    std::vector<float> Pkmid() const;

    /// \returns BamRecord::Pkmean2 with this view's parameters applied
    std::vector<float> Pkmean2() const;

    /// \returns BamRecord::Pkmid2 with this view's parameters applied
    std::vector<float> Pkmid2() const;

    /// \returns BamRecord::PreBaseFrames with this view's parameters applied
    Frames PrebaseFrames() const;

    /// \returns BamRecord::PrePulseFrames with this view's parameters applied
    Frames PrePulseFrames() const;

    /// \returns BamRecord::PulseCalls with this view's parameters applied
    std::string PulseCalls() const;

    /// \returns BamRecord::PulseCallWidth with this view's parameters applied
    Frames PulseCallWidth() const;

    /// \returns BamRecord::PulseWidths with this view's parameters applied
    Frames PulseWidths() const;

    /// \returns BamRecord::Qualities with this view's parameters applied
    QualityValues Qualities() const;

    /// \returns BamRecord::Sequence with this view's parameters applied
    std::string Sequence() const;

    /// \returns BamRecord::StartFrame with this view's parameters applied
    std::vector<uint32_t> StartFrames() const;

    /// \returns BamRecord::SubstitutionQV with this view's parameters applied
    QualityValues SubstitutionQVs() const;

    /// \returns BamRecord::SubstitutionTag with this view's parameters applied
    std::string SubstitutionTags() const;

private:
    const BamRecord& record_;
    Orientation orientation_;
    bool aligned_;
    bool exciseSoftClips_;
    PulseBehavior pulseBehavior_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/BamRecordView.inl"

#endif // BAMRECORDVIEW_H
