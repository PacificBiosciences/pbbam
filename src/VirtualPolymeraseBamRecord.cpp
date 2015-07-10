// Copyright (c) 2015, Pacific Biosciences of California, Inc.
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

// Author: Armin Töpfer

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "pbbam/virtual/VirtualPolymeraseBamRecord.h"
#include "pbbam/virtual/VirtualRegionType.h"
#include "pbbam/virtual/VirtualRegionTypeMap.h"

using namespace PacBio;
using namespace PacBio::BAM;

VirtualPolymeraseBamRecord::VirtualPolymeraseBamRecord(
    std::vector<BamRecord>&& unorderedSources, const BamHeader& header)
    : BamRecord(header)
    , sources_(std::forward<std::vector<BamRecord>>(unorderedSources))
{
    // Sort sources by queryStart
    std::sort(sources_.begin(), sources_.end(),
              [](const BamRecord& l1, const BamRecord& l2)
              { return l1.QueryStart() < l2.QueryStart(); });
    StitchSources();
}

void VirtualPolymeraseBamRecord::StitchSources()
{
    const auto& firstRecord = sources_[0];
    const auto& lastRecord = sources_[sources_.size() - 1];

    // Temporary variables used for stitching
    int accuracy = 0;
    int accuracyCounter = 0;

    std::string   sequence;
    std::string   deletionTag;
    std::string   substitutionTag;
    std::string   labelTag;
    std::string   alternativeLabelTag;
    std::string   pulseCall;

    QualityValues qualities;
    QualityValues deletionQv;
    QualityValues insertionQv;
    QualityValues mergeQv;
    QualityValues substitutionQv;
    QualityValues labelQv;
    QualityValues alternativeLabelQv;

    Frames             ipd;
    Frames             pw;
    Frames             pd;
    Frames             px;
    std::vector<float> pa;
    std::vector<float> pm;

    // Stitch using tmp vars
    for(auto& b : sources_)
    {
        sequence.append(b.Sequence());
        
        MoveAppend(b.Qualities(), qualities);

        if (b.HasReadAccuracy())
        {
            accuracy += b.ReadAccuracy();
            ++accuracyCounter;
        }

        if (b.HasDeletionQV())
            MoveAppend(std::move(b.DeletionQV()), deletionQv);

        if (b.HasInsertionQV())
            MoveAppend(std::move(b.InsertionQV()), insertionQv);

        if (b.HasMergeQV())
            MoveAppend(std::move(b.MergeQV()), mergeQv);

        if (b.HasSubstitutionQV())
            MoveAppend(std::move(b.SubstitutionQV()), substitutionQv);

        if (b.HasLabelQV())
            MoveAppend(std::move(b.LabelQV()), labelQv);

        if (b.HasAltLabelQV())
            MoveAppend(std::move(b.AltLabelQV()), alternativeLabelQv);

        if (b.HasDeletionTag())
            deletionTag.append(std::move(b.DeletionTag()));

        if (b.HasSubstitutionTag())
            substitutionTag.append(std::move(b.SubstitutionTag()));

        if (b.HasLabelTag())
            labelTag.append(std::move(b.LabelTag()));

        if (b.HasAltLabelTag())
            alternativeLabelTag.append(std::move(b.AltLabelTag()));

        if (b.HasPulseCall())
            pulseCall.append(std::move(b.PulseCall()));

        if (b.HasIPD())
            MoveAppend(b.IPDRaw().DataRaw(), ipd.DataRaw());

        if (b.HasPulseWidth())
            MoveAppend(b.PulseWidthRaw().DataRaw(), pw.DataRaw());

        if (b.HasPulseCallWidth())
            MoveAppend(b.PulseCallWidth().DataRaw(), px.DataRaw());

        if (b.HasPrePulseFrames())
            MoveAppend(b.PrePulseFrames().DataRaw(), pd.DataRaw());

        if (b.HasPkmid())
            MoveAppend(b.Pkmid(), pm);

        if (b.HasPkmean())
            MoveAppend(b.Pkmean(), pa);

        if (b.HasScrapType())
        {
            const auto regionType = b.ScrapType();

            if (!HasVirtualRegionType(regionType))
                virtualRegionsMap_[regionType] = std::vector<VirtualRegion>();

            virtualRegionsMap_[regionType].emplace_back(
                regionType, b.QueryStart(), b.QueryEnd());   
        }

        if (b.HasLocalContextFlags())
        {
            std::pair<int, int> barcodes{-1, -1};
            if (b.HasBarcodes())
                barcodes = b.Barcodes();

            constexpr auto regionType = VirtualRegionType::SUBREAD;
            if (!HasVirtualRegionType(regionType))
                virtualRegionsMap_[regionType] = std::vector<VirtualRegion>();

            virtualRegionsMap_[regionType].emplace_back(
                regionType, b.QueryStart(), b.QueryEnd(), b.LocalContextFlags(),
                barcodes.first, barcodes.second);   
        }
    }

    // ReadGroup
    this->ReadGroup(this->header_.ReadGroups()[0]);

    // Avoid division by 0
    if (accuracyCounter > 0)
        this->ReadAccuracy(accuracy / accuracyCounter);

    this->NumPasses(1);

    // All records should contain the same SNR and hole number
    if (firstRecord.HasSignalToNoise())
        this->SignalToNoise(firstRecord.SignalToNoise());
    this->HoleNumber(firstRecord.HoleNumber());

    // QueryStart
    this->QueryStart(firstRecord.QueryStart());
    this->QueryEnd(lastRecord.QueryEnd());
    this->UpdateName();

    std::string qualitiesStr = qualities.Fastq();
    if (sequence.size() == qualitiesStr.size())
        this->Impl().SetSequenceAndQualities(sequence, qualitiesStr);
    else
        this->Impl().SetSequenceAndQualities(sequence);

    // Tags as strings
    if (!deletionTag.empty())
        this->DeletionTag(deletionTag);
    if (!substitutionTag.empty())
        this->SubstitutionTag(substitutionTag);
    if (!labelTag.empty())
        this->LabelTag(labelTag);
    if (!alternativeLabelTag.empty())
        this->AltLabelTag(alternativeLabelTag);
    if (!pulseCall.empty())
        this->PulseCall(pulseCall);

    // QVs
    if (!deletionQv.empty())
        this->DeletionQV(deletionQv);
    if (!insertionQv.empty())
        this->InsertionQV(insertionQv);
    if (!mergeQv.empty())
        this->MergeQV(mergeQv);
    if (!substitutionQv.empty())
        this->SubstitutionQV(substitutionQv);
    if (!labelQv.empty())
        this->LabelQV(labelQv);
    if (!alternativeLabelQv.empty())
        this->AltLabelQV(alternativeLabelQv);

    // 16 bit arrays
    if (!ipd.Data().empty())
        this->IPD(ipd, FrameEncodingType::LOSSLESS);
    if (!pw.Data().empty())
        this->PulseWidth(pw, FrameEncodingType::LOSSLESS);
    if (!pa.empty())
        this->Pkmean(pa);
    if (!pm.empty())
        this->Pkmid(pm);
    if (!pd.Data().empty())
        this->PrePulseFrames(pd, FrameEncodingType::LOSSLESS);
    if (!px.Data().empty())
        this->PulseCallWidth(px, FrameEncodingType::LOSSLESS);

    // Determine HQREGION bases on LQREGIONS
    if (HasVirtualRegionType(VirtualRegionType::LQREGION))
    {
        int beginPos = 0;
        for (const auto& lqregion : virtualRegionsMap_[VirtualRegionType::LQREGION])
        {
            if (lqregion.beginPos - beginPos > 0)
                virtualRegionsMap_[VirtualRegionType::HQREGION].emplace_back(
                    VirtualRegionType::HQREGION, beginPos, lqregion.beginPos);
            beginPos = lqregion.endPos;
        }
    }
}

Frames VirtualPolymeraseBamRecord::IPDV1Frames(Orientation orientation) const
{
    const auto rawFrames = this->IPDRaw(orientation);
    const std::vector<uint8_t> rawData(rawFrames.Data().begin(), rawFrames.Data().end());
    return Frames::Decode(rawData);
}