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

#include "pbbam/Tag.h"
#include <stdexcept>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

Tag::Tag(void)           : data_(),      modifier_(TagModifier::NONE) { }
Tag::Tag(int8_t value)   : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint8_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(int16_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint16_t value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(int32_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint32_t value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(float value)    : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::string& value)      : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<int8_t>& value)   : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<uint8_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<int16_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<uint16_t>& value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<int32_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<uint32_t>& value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<float>& value)    : data_(value), modifier_(TagModifier::NONE) { }

Tag::Tag(int8_t value, const TagModifier mod)
    : data_(value)
    , modifier_(mod)
{
    if (mod == TagModifier::HEX_STRING)
        throw runtime_error("HEX_STRING is not a valid tag modifier for int8_t data. "
                            "It is intended for string-type data only.");
}

Tag::Tag(const std::string& value, const TagModifier mod)
    : data_(value)
    , modifier_(mod)
{
    if (mod == TagModifier::ASCII_CHAR)
        throw runtime_error("ASCII_CHAR is not a valid tag modifier for string-type data. "
                            "To construct an ASCII char tag, use a single-quoted value (e.g. 'X' instead of \"X\")");
}

Tag::Tag(const Tag& other)
    : data_(other.data_)
    , modifier_(other.modifier_)
{ }

Tag::Tag(Tag&& other)
    : data_(std::move(other.data_))
    , modifier_(std::move(other.modifier_))
{ }

Tag::~Tag(void) { }

Tag& Tag::operator=(boost::blank value) { data_ = value; return *this; }
Tag& Tag::operator=(int8_t value)   { data_ = value; return *this; }
Tag& Tag::operator=(uint8_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(int16_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(uint16_t value) { data_ = value; return *this; }
Tag& Tag::operator=(int32_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(uint32_t value) { data_ = value; return *this; }
Tag& Tag::operator=(float value)    { data_ = value; return *this; }
Tag& Tag::operator=(const std::string& value)      { data_ = value; return *this; }
Tag& Tag::operator=(const vector<int8_t>& value)   { data_ = value; return *this; }
Tag& Tag::operator=(const vector<uint8_t>& value)  { data_ = value; return *this; }
Tag& Tag::operator=(const vector<int16_t>& value)  { data_ = value; return *this; }
Tag& Tag::operator=(const vector<uint16_t>& value) { data_ = value; return *this; }
Tag& Tag::operator=(const vector<int32_t>& value)  { data_ = value; return *this; }
Tag& Tag::operator=(const vector<uint32_t>& value) { data_ = value; return *this; }
Tag& Tag::operator=(const vector<float>& value)    { data_ = value; return *this; }

Tag& Tag::operator=(const Tag& other)
{
    data_ = other.data_;
    modifier_ = other.modifier_;
    return *this;
}

Tag& Tag::operator=(Tag&& other)
{
    data_ = std::move(other.data_);
    modifier_ = std::move(other.modifier_);
    return *this;
}
