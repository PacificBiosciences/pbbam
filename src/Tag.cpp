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
/// \file Tag.cpp
/// \brief Defines the Tag class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Tag.h"
#include <cstdint>
#include <stdexcept>

namespace PacBio {
namespace BAM {

Tag::Tag(void)           : data_(),      modifier_(TagModifier::NONE) { }
Tag::Tag(int8_t value)   : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint8_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(int16_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint16_t value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(int32_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint32_t value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(float value)    : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::string& value)      : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::vector<int8_t>& value)   : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::vector<uint8_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::vector<int16_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::vector<uint16_t>& value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::vector<int32_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::vector<uint32_t>& value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::vector<float>& value)    : data_(value), modifier_(TagModifier::NONE) { }

Tag::Tag(int8_t value, const TagModifier mod)
    : data_(value)
    , modifier_(mod)
{
    if (mod == TagModifier::HEX_STRING)
        throw std::runtime_error("HEX_STRING is not a valid tag modifier for int8_t data. "
                                 "It is intended for string-type data only.");
}

Tag::Tag(const std::string& value, const TagModifier mod)
    : data_(value)
    , modifier_(mod)
{
    if (mod == TagModifier::ASCII_CHAR)
        throw std::runtime_error("ASCII_CHAR is not a valid tag modifier for string-type data. "
                                 "To construct an ASCII char tag, use a single-quoted value (e.g. 'X' instead of \"X\")");
}

Tag& Tag::operator=(boost::blank value) { data_ = value; return *this; }
Tag& Tag::operator=(int8_t value)   { data_ = value; return *this; }
Tag& Tag::operator=(uint8_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(int16_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(uint16_t value) { data_ = value; return *this; }
Tag& Tag::operator=(int32_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(uint32_t value) { data_ = value; return *this; }
Tag& Tag::operator=(float value)    { data_ = value; return *this; }
Tag& Tag::operator=(const std::string& value)      { data_ = value; return *this; }
Tag& Tag::operator=(const std::vector<int8_t>& value)   { data_ = value; return *this; }
Tag& Tag::operator=(const std::vector<uint8_t>& value)  { data_ = value; return *this; }
Tag& Tag::operator=(const std::vector<int16_t>& value)  { data_ = value; return *this; }
Tag& Tag::operator=(const std::vector<uint16_t>& value) { data_ = value; return *this; }
Tag& Tag::operator=(const std::vector<int32_t>& value)  { data_ = value; return *this; }
Tag& Tag::operator=(const std::vector<uint32_t>& value) { data_ = value; return *this; }
Tag& Tag::operator=(const std::vector<float>& value)    { data_ = value; return *this; }

} // namespace BAM
} // namespace PacBio
