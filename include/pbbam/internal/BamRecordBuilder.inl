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
/// \file BamRecordBuilder.inl
/// \brief Inline implementations for the BamRecordBuilder class.
//
// Author: Derek Barnett

#include "pbbam/BamRecordBuilder.h"

namespace PacBio {
namespace BAM {

inline BamRecordBuilder& BamRecordBuilder::Bin(const uint32_t bin)
{ core_.bin = bin; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Flag(const uint32_t flag)
{ core_.flag = flag; return *this; }

inline BamRecordBuilder& BamRecordBuilder::InsertSize(const int32_t iSize)
{ core_.isize = iSize; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MapQuality(const uint8_t mapQual)
{ core_.qual = mapQual; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MatePosition(const int32_t pos)
{ core_.mpos = pos; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MateReferenceId(const int32_t id)
{ core_.mtid = id; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Position(const int32_t pos)
{ core_.pos = pos; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Qualities(const std::string& qualities)
{ qualities_ = qualities; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Qualities(std::string&& qualities)
{ qualities_ = std::move(qualities); return *this; }

inline BamRecordBuilder& BamRecordBuilder::ReferenceId(const int32_t id)
{ core_.tid = id; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Tags(const TagCollection& tags)
{ tags_ = tags; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Tags(TagCollection&& tags)
{ tags_ = std::move(tags); return *this; }

} // namespace BAM
} // namespace PacBio
