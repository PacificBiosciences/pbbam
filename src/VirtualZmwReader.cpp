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
/// \file VirtualZmwReader.cpp
/// \brief Implements the VirtualZmwReader class.
//
// Author: Armin Töpfer

#include <stdexcept>

#include "VirtualZmwReader.h"
#include "pbbam/ReadGroupInfo.h"

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

VirtualZmwReader::VirtualZmwReader(const std::string& primaryBamFilepath,
                                   const std::string& scrapsBamFilepath)
    : VirtualZmwReader(primaryBamFilepath, scrapsBamFilepath, PbiFilter{})
{ }

VirtualZmwReader::VirtualZmwReader(const std::string& primaryBamFilepath,
                                   const std::string& scrapsBamFilepath,
                                   const PbiFilter& filter)
{
    primaryBamFile_.reset(new BamFile{ primaryBamFilepath });
    scrapsBamFile_.reset(new BamFile{ scrapsBamFilepath });

    if (filter.IsEmpty()) {
        primaryQuery_.reset(new EntireFileQuery(*primaryBamFile_));
        scrapsQuery_.reset(new EntireFileQuery(*scrapsBamFile_));
    }
    else {
        primaryQuery_.reset(new PbiFilterQuery{ filter, *primaryBamFile_ });
        scrapsQuery_.reset(new PbiFilterQuery{ filter, *scrapsBamFile_ });
    }

    primaryIt_ = (primaryQuery_->begin());
    scrapsIt_ = (scrapsQuery_->begin());

    stitchedHeader_.reset(new BamHeader{ primaryBamFile_->Header().ToSam() });

    // update stitched read group in header
    auto readGroups = stitchedHeader_->ReadGroups();
    if (readGroups.empty())
        throw std::runtime_error("Bam header of the primary bam has no read groups.");
    readGroups[0].ReadType("POLYMERASE");
    readGroups[0].Id(readGroups[0].MovieName(), "POLYMERASE");
    if (readGroups.size() > 1)
    {
        std::vector<ReadGroupInfo> singleGroup;
        singleGroup.emplace_back(std::move(readGroups[0]));
        readGroups = std::move(singleGroup);
        stitchedHeader_->ClearReadGroups();
    }
    stitchedHeader_->ReadGroups(readGroups);
}

VirtualZmwReader::~VirtualZmwReader(void) { }

bool VirtualZmwReader::HasNext(void)
{
    // Return true until both iterators are at the end of the query
    return primaryIt_ != primaryQuery_->end() ||
            scrapsIt_ != scrapsQuery_->end();
}

// This method is not thread safe
VirtualZmwBamRecord VirtualZmwReader::Next(void)
{ return VirtualZmwBamRecord{ NextRaw(), *stitchedHeader_ }; }

std::vector<BamRecord> VirtualZmwReader::NextRaw(void)
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
        currentHoleNumber = std::min((*primaryIt_).HoleNumber(),
                                     (*scrapsIt_).HoleNumber());

    // collect subreads or hqregions
    while (primaryIt_ != primaryQuery_->end() &&
           currentHoleNumber == (*primaryIt_).HoleNumber())
    {
        bamRecordVec.push_back(*primaryIt_++);
    }

    // collect scraps
    while (scrapsIt_ != scrapsQuery_->end() &&
           currentHoleNumber == (*scrapsIt_).HoleNumber())
    {
        bamRecordVec.push_back(*scrapsIt_++);
    }

    return bamRecordVec;
}

BamHeader VirtualZmwReader::PrimaryHeader(void) const
{ return primaryBamFile_->Header(); }

BamHeader VirtualZmwReader::ScrapsHeader(void) const
{ return scrapsBamFile_->Header(); }
