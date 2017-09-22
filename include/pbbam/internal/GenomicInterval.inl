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
/// \file GenomicInterval.inl
/// \brief Inline implementations for the GenomicInterval class.
//
// Author: Derek Barnett

#include "pbbam/GenomicInterval.h"

namespace PacBio {
namespace BAM {

inline std::string GenomicInterval::Name() const
{ return name_; }

inline GenomicInterval& GenomicInterval::Name(const std::string& name)
{ name_ = name; return *this; }

inline PacBio::BAM::Interval<Position> GenomicInterval::Interval() const
{ return interval_; }

inline GenomicInterval& GenomicInterval::Interval(const PacBio::BAM::Interval<Position>& interval)
{ interval_ = interval; return *this; }

inline bool GenomicInterval::IsValid() const
{
    return !name_.empty() &&
           interval_.Start() >= 0 &&
           interval_.Stop()  >= 0 &&
           interval_.IsValid();
}

inline size_t GenomicInterval::Length() const
{ return interval_.Length(); }

inline Position GenomicInterval::Start() const
{ return interval_.Start(); }

inline GenomicInterval& GenomicInterval::Start(const Position start)
{ interval_.Start(start); return *this; }

inline Position GenomicInterval::Stop() const
{ return interval_.Stop(); }

inline GenomicInterval& GenomicInterval::Stop(const Position stop)
{ interval_.Stop(stop); return *this; }

inline bool GenomicInterval::operator==(const GenomicInterval& other) const
{ return name_ == other.name_ && interval_ == other.interval_; }

inline bool GenomicInterval::operator!=(const GenomicInterval& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio
