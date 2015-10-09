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
// Author: Derek Barnett

#ifndef PBIBASICTYPES_H
#define PBIBASICTYPES_H

#include "pbbam/Compare.h"
#include "pbbam/Config.h"
#include <deque>
#include <utility>
#include <vector>

namespace PacBio {
namespace BAM {

//
// Contiguous reads that satisfy a query will be returned as a block.
// This is to help minimize number of seeks (or even unneccesary checks).
//
// An index query can iterate over the lookup result 'IndexResultBlocks' list to
// perform a seek and fetch 'numReads' consecutive records before needing to
// seek again.
//
struct PBBAM_EXPORT IndexResultBlock
{
public:
    IndexResultBlock(void);
    IndexResultBlock(size_t idx, size_t numReads);

public:
    bool operator==(const IndexResultBlock& other) const;
    bool operator!=(const IndexResultBlock& other) const;

public:
    size_t  firstIndex_;
    size_t  numReads_;
    int64_t virtualOffset_;
};

typedef std::deque<IndexResultBlock> IndexResultBlocks;
typedef std::vector<size_t>          IndexList;
typedef std::pair<size_t, size_t>    IndexRange;

inline IndexResultBlock::IndexResultBlock(void)
    : firstIndex_(0)
    , numReads_(0)
    , virtualOffset_(-1)
{ }

inline IndexResultBlock::IndexResultBlock(size_t idx, size_t numReads)
    : firstIndex_(idx)
    , numReads_(numReads)
    , virtualOffset_(-1)
{ }

inline bool IndexResultBlock::operator==(const IndexResultBlock& other) const
{
    return firstIndex_ == other.firstIndex_ &&
           numReads_ == other.numReads_ &&
           virtualOffset_ == other.virtualOffset_;
}

inline bool IndexResultBlock::operator!=(const IndexResultBlock& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio

#endif // PBIBASICTYPES_H
