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
//
// File Description
/// \file VirtualPolymeraseReader.cpp
/// \brief Implements the VirtualPolymeraseReader class.
//
// Author: Armin Töpfer

#include <stdexcept>

#include "pbbam/virtual/VirtualPolymeraseReader.h"
#include "pbbam/ReadGroupInfo.h"

using namespace PacBio;
using namespace PacBio::BAM;

VirtualPolymeraseReader::VirtualPolymeraseReader(
    const std::string& primaryBamFilePath, const std::string& scrapsBamFilePath)
    : primaryBamFilePath_(primaryBamFilePath)
    , scrapsBamFilePath_(scrapsBamFilePath)
{
    primaryBamFile_ = std::unique_ptr<BamFile>(new BamFile(primaryBamFilePath_));
    primaryQuery_   = std::unique_ptr<EntireFileQuery>(new EntireFileQuery(*primaryBamFile_));
    primaryIt_      = primaryQuery_->begin();

    scrapsBamFile_  = std::unique_ptr<BamFile>(new BamFile(scrapsBamFilePath_));
    scrapsQuery_    = std::unique_ptr<EntireFileQuery>(new EntireFileQuery(*scrapsBamFile_));
    scrapsIt_       = scrapsQuery_->begin();

    polyHeader_     = std::unique_ptr<BamHeader>(
                        new BamHeader(primaryBamFile_->Header().ToSam()));

    auto readGroups = polyHeader_->ReadGroups();
    if (readGroups.empty())
        throw std::runtime_error("Bam header of the primary bam has no read groups.");
    readGroups[0].ReadType("POLYMERASE");
    readGroups[0].Id(readGroups[0].MovieName(), "POLYMERASE");
    if (readGroups.size() > 1)
    {
        std::vector<ReadGroupInfo> singleGroup;
        singleGroup.emplace_back(std::move(readGroups[0]));
        readGroups = std::move(singleGroup);
        polyHeader_->ClearReadGroups();
    }
    polyHeader_->ReadGroups(readGroups);
}

bool VirtualPolymeraseReader::HasNext(void)
{
    // Return true until both iterators are at the end of the query
    return primaryIt_ != primaryQuery_->end() || scrapsIt_ != scrapsQuery_->end();
}

// This method is not thread safe
VirtualPolymeraseBamRecord VirtualPolymeraseReader::Next(void)
{
    auto bamRecordVec = NextRaw();
    VirtualPolymeraseBamRecord stitched(std::move(bamRecordVec), *polyHeader_);
    return std::move(stitched);
}

std::vector<BamRecord> VirtualPolymeraseReader::NextRaw(void)
{
    std::vector<BamRecord> bamRecordVec;

    // Current hole number, the smallest of scraps and primary.
    // It can be that the next ZMW is scrap only.
    int currentHoleNumber;
    if (primaryIt_ == primaryQuery_->end())
        currentHoleNumber = (*scrapsIt_).HoleNumber();
    else if (scrapsIt_ == scrapsQuery_->end())
        currentHoleNumber = (*primaryIt_).HoleNumber();
    else
        currentHoleNumber = std::min((*primaryIt_).HoleNumber(), (*scrapsIt_).HoleNumber());

    // collect subreads or hqregions
    while (primaryIt_ != primaryQuery_->end() && currentHoleNumber == (*primaryIt_).HoleNumber())
        bamRecordVec.push_back(*primaryIt_++);

    // collect scraps
    while (scrapsIt_ != scrapsQuery_->end() && currentHoleNumber == (*scrapsIt_).HoleNumber())
        bamRecordVec.push_back(*scrapsIt_++);

    return bamRecordVec;
}

BamHeader VirtualPolymeraseReader::PrimaryHeader(void) const
{ return primaryBamFile_->Header(); }

BamHeader VirtualPolymeraseReader::ScrapsHeader(void) const
{ return scrapsBamFile_->Header(); }
