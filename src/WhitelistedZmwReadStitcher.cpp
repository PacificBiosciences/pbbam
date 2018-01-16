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
/// \file WhitelistedZmwReadStitcher.cpp
/// \brief Implements the WhitelistedZmwReadStitcher class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/virtual/WhitelistedZmwReadStitcher.h"

#include <cassert>
#include <cstdint>

#include "VirtualZmwReader.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiIndexedBamReader.h"

namespace PacBio {
namespace BAM {

struct WhitelistedZmwReadStitcher::WhitelistedZmwReadStitcherPrivate
{
public:
    WhitelistedZmwReadStitcherPrivate(const std::vector<int32_t>& zmwWhitelist,
                                      const std::string& primaryBamFilePath,
                                      const std::string& scrapsBamFilePath)
        : primaryBamFile_(new BamFile{primaryBamFilePath})
        , scrapsBamFile_(new BamFile{scrapsBamFilePath})
        , primaryReader_(new PbiIndexedBamReader{*primaryBamFile_})
        , scrapsReader_(new PbiIndexedBamReader{*scrapsBamFile_})
    {
        // setup new header for stitched data
        polyHeader_ = std::make_unique<BamHeader>(primaryBamFile_->Header().ToSam());
        auto readGroups = polyHeader_->ReadGroups();
        if (readGroups.empty())
            throw std::runtime_error("Bam header of the primary bam has no read groups.");
        readGroups[0].ReadType("POLYMERASE");
        readGroups[0].Id(readGroups[0].MovieName(), "POLYMERASE");
        if (readGroups.size() > 1) {
            std::vector<ReadGroupInfo> singleGroup;
            singleGroup.emplace_back(std::move(readGroups[0]));
            readGroups = std::move(singleGroup);
            polyHeader_->ClearReadGroups();
        }
        polyHeader_->ReadGroups(readGroups);

        // remove ZMWs up front, that are not found in either file
        PreFilterZmws(zmwWhitelist);
    }

    bool HasNext() const { return !zmwWhitelist_.empty(); }

    VirtualZmwBamRecord Next()
    {
        auto bamRecordVec = NextRaw();
        VirtualZmwBamRecord stitched(std::move(bamRecordVec), *polyHeader_);
        return stitched;
    }

    std::vector<BamRecord> NextRaw()
    {
        auto result = std::vector<BamRecord>{};
        if (!HasNext()) return result;

        const auto& zmw = zmwWhitelist_.front();
        primaryReader_->Filter(PbiZmwFilter{zmw});
        scrapsReader_->Filter(PbiZmwFilter{zmw});

        auto record = BamRecord{};
        while (primaryReader_->GetNext(record))
            result.push_back(record);
        while (scrapsReader_->GetNext(record))
            result.push_back(record);

        zmwWhitelist_.pop_front();
        return result;
    }

    BamHeader PrimaryHeader() const { return primaryBamFile_->Header(); }

    BamHeader ScrapsHeader() const { return scrapsBamFile_->Header(); }

private:
    std::unique_ptr<BamFile> primaryBamFile_;
    std::unique_ptr<BamFile> scrapsBamFile_;
    std::unique_ptr<PbiIndexedBamReader> primaryReader_;
    std::unique_ptr<PbiIndexedBamReader> scrapsReader_;
    std::unique_ptr<BamHeader> polyHeader_;
    std::deque<int32_t> zmwWhitelist_;

private:
    void PreFilterZmws(const std::vector<int32_t>& zmwWhitelist)
    {
        // fetch input ZMWs
        const PbiRawData primaryIndex(primaryBamFile_->PacBioIndexFilename());
        const PbiRawData scrapsIndex(scrapsBamFile_->PacBioIndexFilename());
        const auto& primaryZmws = primaryIndex.BasicData().holeNumber_;
        const auto& scrapsZmws = scrapsIndex.BasicData().holeNumber_;

        // toss them all into a set (for uniqueness & lookup here soon)
        std::set<int32_t> inputZmws;
        for (const auto& zmw : primaryZmws)
            inputZmws.insert(zmw);
        for (const auto& zmw : scrapsZmws)
            inputZmws.insert(zmw);

        // check our requested whitelist against files' ZMWs, keep if found
        const auto inputEnd = inputZmws.cend();
        for (const int32_t zmw : zmwWhitelist) {
            if (inputZmws.find(zmw) != inputEnd) zmwWhitelist_.push_back(zmw);
        }
    }
};

// --------------------------------
// ZmwReadStitcher implementation
// --------------------------------

WhitelistedZmwReadStitcher::WhitelistedZmwReadStitcher(const std::vector<int32_t>& zmwWhitelist,
                                                       const std::string& primaryBamFilePath,
                                                       const std::string& scrapsBamFilePath)
    : d_(new WhitelistedZmwReadStitcherPrivate(zmwWhitelist, primaryBamFilePath, scrapsBamFilePath))
{
}

WhitelistedZmwReadStitcher::~WhitelistedZmwReadStitcher() {}

bool WhitelistedZmwReadStitcher::HasNext() const { return d_->HasNext(); }

VirtualZmwBamRecord WhitelistedZmwReadStitcher::Next() { return d_->Next(); }

std::vector<BamRecord> WhitelistedZmwReadStitcher::NextRaw() { return d_->NextRaw(); }

BamHeader WhitelistedZmwReadStitcher::PrimaryHeader() const { return d_->PrimaryHeader(); }

BamHeader WhitelistedZmwReadStitcher::ScrapsHeader() const { return d_->ScrapsHeader(); }

}  // namespace BAM
}  // namespace PacBio
