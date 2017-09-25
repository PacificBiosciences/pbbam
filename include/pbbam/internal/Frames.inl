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
/// \file Frames.inl
/// \brief Inline implementations for the Frames class.
//
// Author: Derek Barnett

#include "pbbam/Frames.h"

namespace PacBio {
namespace BAM {

inline const std::vector<uint16_t>& Frames::Data() const
{ return data_; }

inline std::vector<uint16_t>& Frames::DataRaw()
{ return data_; }

inline std::vector<uint8_t> Frames::Encode() const
{ return Frames::Encode(data_); }

inline Frames& Frames::Data(const std::vector<uint16_t>& frames)
{ data_ = frames; return *this; }

inline Frames& Frames::Data(std::vector<uint16_t>&& frames)
{ data_ = std::move(frames); return *this; }

inline std::vector<uint16_t>::const_iterator Frames::begin() const
{ return data_.begin(); }

inline std::vector<uint16_t>::iterator Frames::begin()
{ return data_.begin(); }

inline std::vector<uint16_t>::const_iterator Frames::cbegin() const
{ return data_.cbegin(); }

inline std::vector<uint16_t>::const_iterator Frames::cend() const
{ return data_.cend(); }

inline std::vector<uint16_t>::const_iterator Frames::end() const
{ return data_.end(); }

inline std::vector<uint16_t>::iterator Frames::end()
{ return data_.end(); }

inline size_t Frames::size() const
{ return data_.size(); }

inline bool Frames::empty() const
{ return data_.empty(); }

inline bool Frames::operator==(const Frames& other) const
{ return data_ == other.data_; }

inline bool Frames::operator!=(const Frames& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio
