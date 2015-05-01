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

#include "pbbam/Tag.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

struct AsciiConvertVisitor : public boost::static_visitor<char>
{
    char operator() (const boost::blank&) const { return '\0'; }
    char operator() (const int8_t& x) const     { return ( x >=33 && x <= 127 ? static_cast<char>(x) : '\0' ); }
    char operator() (const uint8_t& x) const    { return ( x >=33 && x <= 127 ? static_cast<char>(x) : '\0' ); }
    char operator() (const int16_t& x) const    { return ( x >=33 && x <= 127 ? static_cast<char>(x) : '\0' ); }
    char operator() (const uint16_t& x) const   { return ( x >=33 && x <= 127 ? static_cast<char>(x) : '\0' ); }
    char operator() (const int32_t& x) const    { return ( x >=33 && x <= 127 ? static_cast<char>(x) : '\0' ); }
    char operator() (const uint32_t& x) const   { return ( x >=33 && x <= 127 ? static_cast<char>(x) : '\0' ); }
    char operator() (const float&) const            { return '\0'; }
    char operator() (const string&) const           { return '\0'; }
    char operator() (const vector<int8_t>&) const   { return '\0'; }
    char operator() (const vector<uint8_t>&) const  { return '\0'; }
    char operator() (const vector<int16_t>&) const  { return '\0'; }
    char operator() (const vector<uint16_t>&) const { return '\0'; }
    char operator() (const vector<int32_t>&) const  { return '\0'; }
    char operator() (const vector<uint32_t>&) const { return '\0'; }
    char operator() (const vector<float>&) const    { return '\0'; }
};

struct IsEqualVisitor : public boost::static_visitor<bool>
{
    template <typename T, typename U>
    bool operator() (const T&, const U&) const
    {
        // maybe allow conversions down the road?
        // but for now, just fail if types are different
        return false;
    }

    template <typename T>
    bool operator() (const T& lhs, const T& rhs) const
    { return IsEqual(lhs, rhs); }

private:
    bool IsEqual(const boost::blank&, const boost::blank&) const { return true; }
    bool IsEqual(const int8_t& lhs, const int8_t& rhs) const     { return lhs == rhs; }
    bool IsEqual(const uint8_t& lhs, const uint8_t& rhs) const   { return lhs == rhs; }
    bool IsEqual(const int16_t& lhs, const int16_t& rhs) const   { return lhs == rhs; }
    bool IsEqual(const uint16_t& lhs, const uint16_t& rhs) const { return lhs == rhs; }
    bool IsEqual(const int32_t& lhs, const int32_t& rhs) const   { return lhs == rhs; }
    bool IsEqual(const uint32_t& lhs, const uint32_t& rhs) const { return lhs == rhs; }
    bool IsEqual(const float& lhs, const float& rhs) const       { return lhs == rhs; }
    bool IsEqual(const string& lhs, const string& rhs) const     { return lhs == rhs; }
    bool IsEqual(const vector<int8_t>& lhs, const vector<int8_t>& rhs) const     { return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin()); }
    bool IsEqual(const vector<uint8_t>& lhs, const vector<uint8_t>& rhs) const   { return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin()); }
    bool IsEqual(const vector<int16_t>& lhs, const vector<int16_t>& rhs) const   { return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin()); }
    bool IsEqual(const vector<uint16_t>& lhs, const vector<uint16_t>& rhs) const { return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin()); }
    bool IsEqual(const vector<int32_t>& lhs, const vector<int32_t>& rhs) const   { return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin()); }
    bool IsEqual(const vector<uint32_t>& lhs, const vector<uint32_t>& rhs) const { return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin()); }
    bool IsEqual(const vector<float>& lhs, const vector<float>& rhs) const       { return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin()); }
};

struct TypenameVisitor : public boost::static_visitor<string>
{
    string operator() (const boost::blank&) const     { return "none"; }
    string operator() (const int8_t&) const           { return "int8_t"; }
    string operator() (const uint8_t&) const          { return "uint8_t"; }
    string operator() (const int16_t&) const          { return "int16_t"; }
    string operator() (const uint16_t&) const         { return "uint16_t"; }
    string operator() (const int32_t&) const          { return "int32_t"; }
    string operator() (const uint32_t&) const         { return "uint32_t"; }
    string operator() (const float&) const            { return "float"; }
    string operator() (const string&) const           { return "string"; }
    string operator() (const vector<int8_t>&) const   { return "vector<int8_t>"; }
    string operator() (const vector<uint8_t>&) const  { return "vector<uint8_t>"; }
    string operator() (const vector<int16_t>&) const  { return "vector<int16_t>"; }
    string operator() (const vector<uint16_t>&) const { return "vector<uint16_t>"; }
    string operator() (const vector<int32_t>&) const  { return "vector<int32_t>"; }
    string operator() (const vector<uint32_t>&) const { return "vector<uint32_t>"; }
    string operator() (const vector<float>&) const    { return "vector<float>"; }
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

// ------------------------
// Tag implementation
// ------------------------

Tag::Tag(void) : data_(), modifier_(TagModifier::NONE) { }
Tag::Tag(int8_t value)   : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint8_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(int16_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint16_t value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(int32_t value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(uint32_t value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(float value)    : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const std::string& value)           : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<int8_t>& value)   : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<uint8_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<int16_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<uint16_t>& value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<int32_t>& value)  : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<uint32_t>& value) : data_(value), modifier_(TagModifier::NONE) { }
Tag::Tag(const vector<float>& value)    : data_(value), modifier_(TagModifier::NONE) { }

Tag::Tag(const Tag& other)
    : data_(other.data_)
    , modifier_(other.modifier_)
{ }

Tag::Tag(Tag&& other)
    : data_(std::move(other.data_))
    , modifier_(std::move(other.modifier_))
{}

Tag::~Tag(void) { }

Tag& Tag::operator=(boost::blank value) { data_ = value; return *this; }
Tag& Tag::operator=(int8_t value)   { data_ = value; return *this; }
Tag& Tag::operator=(uint8_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(int16_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(uint16_t value) { data_ = value; return *this; }
Tag& Tag::operator=(int32_t value)  { data_ = value; return *this; }
Tag& Tag::operator=(uint32_t value) { data_ = value; return *this; }
Tag& Tag::operator=(float value)    { data_ = value; return *this; }
Tag& Tag::operator=(const std::string& value)           { data_ = value; return *this; }
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

bool Tag::operator== (const Tag& other) const
{
    return boost::apply_visitor(internal::IsEqualVisitor(), data_, other.data_) &&
           (modifier_ == other.modifier_) ;
}

bool Tag::operator!= (const Tag& other) const
{ return !(*this == other); }

bool Tag::HasModifier(const TagModifier m) const
{
    // we just allow one at a time (for now at least)
    return modifier_ == m;
}

TagModifier Tag::Modifier(void) const
{ return modifier_; }

Tag& Tag::Modifier(const TagModifier m)
{ modifier_ = m; return *this; }

char Tag::ToAscii(void) const
{ return boost::apply_visitor(internal::AsciiConvertVisitor(), data_); }

TagDataType Tag::Type(void) const
{ return TagDataType( data_.which() ); }

string Tag::Typename(void) const
{ return boost::apply_visitor(internal::TypenameVisitor(), data_); }
