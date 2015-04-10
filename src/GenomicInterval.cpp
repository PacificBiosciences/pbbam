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

#include "pbbam/GenomicInterval.h"
#include "AssertUtils.h"
#include <htslib/hts.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

GenomicInterval::GenomicInterval(void)
    : id_(-1)
{ }

GenomicInterval::GenomicInterval(const int id,
                                 const Position& start,
                                 const Position& stop)
    : id_(id)
    , interval_(start, stop)
{ }

GenomicInterval::GenomicInterval(const GenomicInterval& other)
    : id_(other.id_)
    , interval_(other.interval_)
{ }

bool GenomicInterval::CoveredBy(const GenomicInterval& other) const
{
    if (id_ != other.id_)
        return false;
    return interval_.CoveredBy(other.interval_);
}

bool GenomicInterval::Covers(const GenomicInterval& other) const
{
    if (id_ != other.id_)
        return false;
    return interval_.Covers(other.interval_);
}

bool GenomicInterval::Intersects(const GenomicInterval& other) const
{
    if (id_ != other.id_)
        return false;
    return interval_.Intersects(other.interval_);
}

size_t GenomicInterval::Length(void) const
{
    return interval_.Length();
}
