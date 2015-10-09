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

#include "PbiIndexedBamReader.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

PbiIndexedBamReader::PbiIndexedBamReader(const PbiFilter& filter,
                                         const std::string& filename)
    : PbiIndexedBamReader(filter, BamFile(filename))
{ }

PbiIndexedBamReader::PbiIndexedBamReader(const PbiFilter& filter,
                                         const BamFile& bamFile)
    : BamReader(bamFile)
    , index_(bamFile_.PacBioIndexFilename())
    , currentBlockReadCount_(0)
{
    Filter(filter);
}

PbiIndexedBamReader::PbiIndexedBamReader(const PbiFilter& filter,
                                         BamFile&& bamFile)
    : BamReader(std::move(bamFile))
    , index_(bamFile_.PacBioIndexFilename())
    , currentBlockReadCount_(0)
{
    Filter(filter);
}

int PbiIndexedBamReader::ReadRawData(BGZF* bgzf, bam1_t* b)
{
    // no data to fetch, return false
    if (blocks_.empty())
        return -1; // "EOF"

    // if on new block, seek to its first record
    if (currentBlockReadCount_ == 0) {
        auto seekResult = bgzf_seek(bgzf, blocks_.at(0).virtualOffset_, SEEK_SET);
        if (seekResult == -1)
            throw std::runtime_error("could not seek in BAM file");
    }

    // read next record
    auto result = bam_read1(bgzf, b);

    // update counters. if block finished, pop & reset
    ++currentBlockReadCount_;
    if (currentBlockReadCount_ == blocks_.at(0).numReads_) {
        blocks_.pop_front();
        currentBlockReadCount_ = 0;
    }

    return result;
}

const PbiFilter& PbiIndexedBamReader::Filter(void) const
{ return filter_; }

PbiIndexedBamReader& PbiIndexedBamReader::Filter(const PbiFilter& filter)
{
    // store request & reset counters
    filter_ = filter;
    currentBlockReadCount_ = 0;
    blocks_.clear();

    // query for index blocks
    auto indexList = filter.Lookup(index_);
    blocks_ = mergedIndexBlocks(indexList);
    index_.BasicData().ApplyOffsets(blocks_);
    return *this;
}
