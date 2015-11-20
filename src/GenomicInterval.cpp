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
/// \file GenomicInterval.cpp
/// \brief Implements the GenomicInterval class.
//
// Author: Derek Barnett

#include "pbbam/GenomicInterval.h"
#include "AssertUtils.h"
#include "StringUtils.h"
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <stdexcept>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// returns sequence name & sets begin/end, from input regionString
string parseRegionString(const string& reg,
                         PacBio::BAM::Position* begin,
                         PacBio::BAM::Position* end)
{
    const vector<string> parts = internal::Split(reg, ':');
    if (parts.empty() || parts.size() > 2)
        throw std::runtime_error("malformed region string");

    // given name only, default min,max intervals
    if (parts.size() == 1) {
        *begin = 0;
        *end = 1<<29;
    }

    // parse interval from input
    else if (parts.size() == 2) {
        const vector<string> intervalParts = internal::Split(parts.at(1), '-');
        if (intervalParts.empty() || intervalParts.size() >2 )
            throw std::runtime_error("malformed region string");
        *begin = std::stoi(intervalParts.at(0));
        *end   = std::stoi(intervalParts.at(1));
    }

    return parts.at(0);
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

GenomicInterval::GenomicInterval(void) { }

GenomicInterval::GenomicInterval(const std::string& name,
                                 const Position& start,
                                 const Position& stop)
    : name_(name)
    , interval_(start, stop)
{ }

GenomicInterval::GenomicInterval(const string& samtoolsRegionString)
{
    Position begin;
    Position end;
    name_ = internal::parseRegionString(samtoolsRegionString, &begin, &end);
    interval_ = PacBio::BAM::Interval<Position>(begin, end);
}

GenomicInterval::GenomicInterval(const GenomicInterval& other)
    : name_(other.name_)
    , interval_(other.interval_)
{ }

GenomicInterval& GenomicInterval::operator=(const GenomicInterval& other)
{
    name_ = other.name_;
    interval_ = other.interval_;
    return *this;
}

bool GenomicInterval::CoveredBy(const GenomicInterval& other) const
{
    if (name_ != other.name_)
        return false;
    return interval_.CoveredBy(other.interval_);
}

bool GenomicInterval::Covers(const GenomicInterval& other) const
{
    if (name_ != other.name_)
        return false;
    return interval_.Covers(other.interval_);
}

bool GenomicInterval::Intersects(const GenomicInterval& other) const
{
    if (name_ != other.name_)
        return false;
    return interval_.Intersects(other.interval_);
}
