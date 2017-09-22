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
/// \file Tag.inl
/// \brief Inline implementations for the Tag class.
//
// Author: Derek Barnett

#include "pbbam/Tag.h"
#include <boost/numeric/conversion/cast.hpp>
#include <iostream>

namespace PacBio {
namespace BAM {
namespace internal {

template<typename T>
inline bool InAsciiRange(const T& x)
{ return (x >=33 && x <= 127); }

struct AsciiConvertVisitor : public boost::static_visitor<char>
{
    // only valid for numeric types - maybe even more restrictive?
    char operator() (const int8_t& x) const   { return Helper(x); }
    char operator() (const uint8_t& x) const  { return Helper(x); }
    char operator() (const int16_t& x) const  { return Helper(x); }
    char operator() (const uint16_t& x) const { return Helper(x); }
    char operator() (const int32_t& x) const  { return Helper(x); }
    char operator() (const uint32_t& x) const { return Helper(x); }

    // anything else always throws
    template<typename T>
    char operator()(const T&) const
    { throw std::runtime_error("conversion not supported"); return 0; }

private:
    template<typename T>
    char Helper(const T& x) const
    {
        if (!InAsciiRange(x))
            throw std::runtime_error("not valid ASCII");
        return static_cast<char>(x);
    }
};

template<typename DesiredType>
struct NumericConvertVisitor : public boost::static_visitor<DesiredType>
{
    // only valid for integral types
    DesiredType operator() (const int8_t& x) const   { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator() (const uint8_t& x) const  { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator() (const int16_t& x) const  { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator() (const uint16_t& x) const { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator() (const int32_t& x) const  { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator() (const uint32_t& x) const { return boost::numeric_cast<DesiredType>(x); }

    // anything else always throws
    template<typename T> DesiredType operator()(const T& t) const
    {
        const std::string from = typeid(t).name();
        const std::string to   = typeid(DesiredType).name();
        const std::string msg  = std::string("conversion not supported: ") + from + " -> " + to;
        throw std::runtime_error(msg);
        return 0;
    }
};

using ToInt8ConvertVisitor   = NumericConvertVisitor<int8_t>;
using ToUInt8ConvertVisitor  = NumericConvertVisitor<uint8_t>;
using ToInt16ConvertVisitor  = NumericConvertVisitor<int16_t>;
using ToUInt16ConvertVisitor = NumericConvertVisitor<uint16_t>;
using ToInt32ConvertVisitor  = NumericConvertVisitor<int32_t>;
using ToUInt32ConvertVisitor = NumericConvertVisitor<uint32_t>;

struct IsEqualVisitor : public boost::static_visitor<bool>
{
    template <typename T, typename U>
    bool operator() (const T&, const U&) const
    {
        // maybe allow conversions down the road?
        // but for now, just fail if types are different
        return false;
    }

    bool operator() (const boost::blank&, const boost::blank&) const
    { return true; }

    template <typename T>
    bool operator() (const T& lhs, const T& rhs) const
    { return lhs == rhs; }
};

struct TypenameVisitor : public boost::static_visitor<std::string>
{
    std::string operator() (const boost::blank&) const     { return "none"; }
    std::string operator() (const int8_t&) const           { return "int8_t"; }
    std::string operator() (const uint8_t&) const          { return "uint8_t"; }
    std::string operator() (const int16_t&) const          { return "int16_t"; }
    std::string operator() (const uint16_t&) const         { return "uint16_t"; }
    std::string operator() (const int32_t&) const          { return "int32_t"; }
    std::string operator() (const uint32_t&) const         { return "uint32_t"; }
    std::string operator() (const float&) const            { return "float"; }
    std::string operator() (const std::string&) const      { return "string"; }
    std::string operator() (const std::vector<int8_t>&) const   { return "vector<int8_t>"; }
    std::string operator() (const std::vector<uint8_t>&) const  { return "vector<uint8_t>"; }
    std::string operator() (const std::vector<int16_t>&) const  { return "vector<int16_t>"; }
    std::string operator() (const std::vector<uint16_t>&) const { return "vector<uint16_t>"; }
    std::string operator() (const std::vector<int32_t>&) const  { return "vector<int32_t>"; }
    std::string operator() (const std::vector<uint32_t>&) const { return "vector<uint32_t>"; }
    std::string operator() (const std::vector<float>&) const    { return "vector<float>"; }
};

} // namespace internal

inline bool Tag::operator== (const Tag& other) const
{
    return boost::apply_visitor(internal::IsEqualVisitor(), data_, other.data_) &&
           (modifier_ == other.modifier_) ;
}

inline bool Tag::operator!= (const Tag& other) const
{ return !(*this == other); }

inline bool Tag::HasModifier(const TagModifier m) const
{
    // we just allow one at a time (for now at least)
    return modifier_ == m;
}

inline bool Tag::IsNull() const
{ return Type() == TagDataType::INVALID; }

inline bool Tag::IsInt8() const
{ return Type() == TagDataType::INT8; }

inline bool Tag::IsUInt8() const
{ return Type() == TagDataType::UINT8; }

inline bool Tag::IsInt16() const
{ return Type() == TagDataType::INT16; }

inline bool Tag::IsUInt16() const
{ return Type() == TagDataType::UINT16; }

inline bool Tag::IsInt32() const
{ return Type() == TagDataType::INT32; }

inline bool Tag::IsUInt32() const
{ return Type() == TagDataType::UINT32; }

inline bool Tag::IsFloat() const
{ return Type() == TagDataType::FLOAT; }

inline bool Tag::IsString() const
{ return Type() == TagDataType::STRING; }

inline bool Tag::IsHexString() const
{ return IsString() && modifier_ == TagModifier::HEX_STRING; }

inline bool Tag::IsInt8Array() const
{ return Type() == TagDataType::INT8_ARRAY; }

inline bool Tag::IsUInt8Array() const
{ return Type() == TagDataType::UINT8_ARRAY; }

inline bool Tag::IsInt16Array() const
{ return Type() == TagDataType::INT16_ARRAY; }

inline bool Tag::IsUInt16Array() const
{ return Type() == TagDataType::UINT16_ARRAY; }

inline bool Tag::IsInt32Array() const
{ return Type() == TagDataType::INT32_ARRAY; }

inline bool Tag::IsUInt32Array() const
{ return Type() == TagDataType::UINT32_ARRAY; }

inline bool Tag::IsFloatArray() const
{ return Type() == TagDataType::FLOAT_ARRAY; }

inline bool Tag::IsSignedInt() const
{ return IsInt8() || IsInt16() || IsInt32(); }

inline bool Tag::IsUnsignedInt() const
{ return IsUInt8() || IsUInt16() || IsUInt32(); }

inline bool Tag::IsIntegral() const
{ return IsSignedInt() || IsUnsignedInt(); }

inline bool Tag::IsNumeric() const
{ return IsIntegral() || IsFloat(); }

inline bool Tag::IsSignedArray() const
{ return IsInt8Array() || IsInt16Array() || IsInt32Array(); }

inline bool Tag::IsUnsignedArray() const
{ return IsUInt8Array() || IsUInt16Array() || IsUInt32Array(); }

inline bool Tag::IsIntegralArray() const
{ return IsSignedArray() || IsUnsignedArray(); }

inline bool Tag::IsArray() const
{ return IsIntegralArray() || IsFloatArray(); }

inline TagModifier Tag::Modifier() const
{ return modifier_; }

inline Tag& Tag::Modifier(const TagModifier m)
{ modifier_ = m; return *this; }

inline char Tag::ToAscii() const
{ return boost::apply_visitor(internal::AsciiConvertVisitor(), data_); }

inline int8_t Tag::ToInt8() const
{
    if (IsInt8())
        return boost::get<int8_t>(data_);
    return boost::apply_visitor(internal::ToInt8ConvertVisitor(), data_);
}

inline uint8_t Tag::ToUInt8() const
{
    if (IsUInt8())
        return boost::get<uint8_t>(data_);
    return boost::apply_visitor(internal::ToUInt8ConvertVisitor(), data_);
}

inline int16_t Tag::ToInt16() const
{
    if (IsInt16())
        return boost::get<int16_t>(data_);
    return boost::apply_visitor(internal::ToInt16ConvertVisitor(), data_);
}

inline uint16_t Tag::ToUInt16() const
{
    if (IsUInt16())
        return boost::get<uint16_t>(data_);
    return boost::apply_visitor(internal::ToUInt16ConvertVisitor(), data_);
}

inline int32_t Tag::ToInt32() const
{
    if (IsInt32())
        return boost::get<int32_t>(data_);
    return boost::apply_visitor(internal::ToInt32ConvertVisitor(), data_);
}

inline uint32_t Tag::ToUInt32() const
{
    if (IsUInt32())
        return boost::get<uint32_t>(data_);
    return boost::apply_visitor(internal::ToUInt32ConvertVisitor(), data_);
}

inline float Tag::ToFloat() const
{ return boost::get<float>(data_); }

inline std::string Tag::ToString() const
{ return boost::get<std::string>(data_); }

inline std::vector<int8_t> Tag::ToInt8Array() const
{ return boost::get< std::vector<int8_t> >(data_); }

inline std::vector<uint8_t> Tag::ToUInt8Array() const
{ return boost::get< std::vector<uint8_t> >(data_); }

inline std::vector<int16_t> Tag::ToInt16Array() const
{ return boost::get< std::vector<int16_t> >(data_); }

inline std::vector<uint16_t> Tag::ToUInt16Array() const
{ return boost::get< std::vector<uint16_t> >(data_); }

inline std::vector<int32_t> Tag::ToInt32Array() const
{ return boost::get< std::vector<int32_t> >(data_); }

inline std::vector<uint32_t> Tag::ToUInt32Array() const
{ return boost::get< std::vector<uint32_t> >(data_); }

inline std::vector<float> Tag::ToFloatArray() const
{ return boost::get< std::vector<float> >(data_); }

inline TagDataType Tag::Type() const
{ return TagDataType(data_.which()  ); }

inline std::string Tag::Typename() const
{ return boost::apply_visitor(internal::TypenameVisitor(), data_); }

} // namespace BAM
} // namespace PacBio
