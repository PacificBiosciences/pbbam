// File Description
/// \file Tag.cpp
/// \brief Implements the Tag class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Tag.h"

#include <cassert>

#include <iostream>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

namespace PacBio {
namespace BAM {
namespace {

template <typename T>
bool InAsciiRange(const T x)
{
    return (x >= 33 && x <= 127);
}

struct AsciiConvertVisitor : public boost::static_visitor<char>
{
    // only valid for numeric types - maybe even more restrictive?
    char operator()(const int8_t& x) const { return Helper(x); }
    char operator()(const uint8_t& x) const { return Helper(x); }
    char operator()(const int16_t& x) const { return Helper(x); }
    char operator()(const uint16_t& x) const { return Helper(x); }
    char operator()(const int32_t& x) const { return Helper(x); }
    char operator()(const uint32_t& x) const { return Helper(x); }

    // anything else always throws
    template <typename T>
    char operator()(const T&) const
    {
        throw std::runtime_error{"[pbbam] tag ERROR: cannot convert to ASCII"};
        return 0;
    }

private:
    template <typename T>
    char Helper(const T& x) const
    {
        if (!InAsciiRange(x))
            throw std::runtime_error{"[pbbam] tag ERROR: char is outside valid ASCII range"};
        return static_cast<char>(x);
    }
};

template <typename DesiredType>
struct NumericConvertVisitor : public boost::static_visitor<DesiredType>
{
    // only valid for integral types
    DesiredType operator()(const int8_t& x) const { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator()(const uint8_t& x) const { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator()(const int16_t& x) const { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator()(const uint16_t& x) const { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator()(const int32_t& x) const { return boost::numeric_cast<DesiredType>(x); }
    DesiredType operator()(const uint32_t& x) const { return boost::numeric_cast<DesiredType>(x); }

    // anything else always throws
    template <typename T>
    DesiredType operator()(const T& t) const
    {
        const std::string from = typeid(t).name();
        const std::string to = typeid(DesiredType).name();
        const std::string msg = "[pbbam] tag ERROR: cannot convert type " + from + " to " + to;
        throw std::runtime_error(msg);
        return 0;
    }
};

using ToInt8ConvertVisitor = NumericConvertVisitor<int8_t>;
using ToUInt8ConvertVisitor = NumericConvertVisitor<uint8_t>;
using ToInt16ConvertVisitor = NumericConvertVisitor<int16_t>;
using ToUInt16ConvertVisitor = NumericConvertVisitor<uint16_t>;
using ToInt32ConvertVisitor = NumericConvertVisitor<int32_t>;
using ToUInt32ConvertVisitor = NumericConvertVisitor<uint32_t>;

struct IsEqualVisitor : public boost::static_visitor<bool>
{
    template <typename T, typename U>
    bool operator()(const T&, const U&) const
    {
        // maybe allow conversions down the road?
        // but for now, just fail if types are different
        return false;
    }

    bool operator()(const boost::blank&, const boost::blank&) const { return true; }

    template <typename T>
    bool operator()(const T& lhs, const T& rhs) const
    {
        return lhs == rhs;
    }
};

struct TypenameVisitor : public boost::static_visitor<std::string>
{
    std::string operator()(const boost::blank&) const { return "none"; }
    std::string operator()(const int8_t&) const { return "int8_t"; }
    std::string operator()(const uint8_t&) const { return "uint8_t"; }
    std::string operator()(const int16_t&) const { return "int16_t"; }
    std::string operator()(const uint16_t&) const { return "uint16_t"; }
    std::string operator()(const int32_t&) const { return "int32_t"; }
    std::string operator()(const uint32_t&) const { return "uint32_t"; }
    std::string operator()(const float&) const { return "float"; }
    std::string operator()(const std::string&) const { return "string"; }
    std::string operator()(const std::vector<int8_t>&) const { return "vector<int8_t>"; }
    std::string operator()(const std::vector<uint8_t>&) const { return "vector<uint8_t>"; }
    std::string operator()(const std::vector<int16_t>&) const { return "vector<int16_t>"; }
    std::string operator()(const std::vector<uint16_t>&) const { return "vector<uint16_t>"; }
    std::string operator()(const std::vector<int32_t>&) const { return "vector<int32_t>"; }
    std::string operator()(const std::vector<uint32_t>&) const { return "vector<uint32_t>"; }
    std::string operator()(const std::vector<float>&) const { return "vector<float>"; }
};

}  // namespace

static_assert(std::is_copy_constructible<Tag>::value, "Tag(const Tag&) is not = default");
static_assert(std::is_copy_assignable<Tag>::value, "Tag& operator=(const Tag&) is not = default");

static_assert(std::is_nothrow_move_constructible<Tag>::value, "Tag(Tag&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<Tag>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

Tag::Tag(int8_t value) : data_{value} {}
Tag::Tag(uint8_t value) : data_{value} {}
Tag::Tag(int16_t value) : data_{value} {}
Tag::Tag(uint16_t value) : data_{value} {}
Tag::Tag(int32_t value) : data_{value} {}
Tag::Tag(uint32_t value) : data_{value} {}
Tag::Tag(float value) : data_{value} {}
Tag::Tag(std::string value) : data_{std::move(value)} {}
Tag::Tag(std::vector<int8_t> value) : data_{std::move(value)} {}
Tag::Tag(std::vector<uint8_t> value) : data_{std::move(value)} {}
Tag::Tag(std::vector<int16_t> value) : data_{std::move(value)} {}
Tag::Tag(std::vector<uint16_t> value) : data_{std::move(value)} {}
Tag::Tag(std::vector<int32_t> value) : data_{std::move(value)} {}
Tag::Tag(std::vector<uint32_t> value) : data_{std::move(value)} {}
Tag::Tag(std::vector<float> value) : data_{std::move(value)} {}

Tag::Tag(int8_t value, const TagModifier mod) : data_{value}, modifier_(mod)
{
    if (mod == TagModifier::HEX_STRING)
        throw std::runtime_error{
            "[pbbam] tag ERROR: HEX_STRING is not a valid tag modifier for int8_t data. "
            "It is intended for string-type data only."};
}

Tag::Tag(std::string value, TagModifier mod) : data_{std::move(value)}, modifier_{mod}
{
    if (mod == TagModifier::ASCII_CHAR)
        throw std::runtime_error{
            "[pbbam] tag ERROR: ASCII_CHAR is not a valid tag modifier for string-type data. "
            "To construct an ASCII char tag, use a single-quoted value (e.g. 'X' instead of "
            "\"X\")"};
}

Tag& Tag::operator=(boost::blank value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(int8_t value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(uint8_t value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(int16_t value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(uint16_t value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(int32_t value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(uint32_t value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(float value)
{
    data_ = value;
    return *this;
}

Tag& Tag::operator=(std::string value)
{
    data_ = std::move(value);
    return *this;
}

Tag& Tag::operator=(std::vector<int8_t> value)
{
    data_ = std::move(value);
    return *this;
}

Tag& Tag::operator=(std::vector<uint8_t> value)
{
    data_ = std::move(value);
    return *this;
}

Tag& Tag::operator=(std::vector<int16_t> value)
{
    data_ = std::move(value);
    return *this;
}

Tag& Tag::operator=(std::vector<uint16_t> value)
{
    data_ = std::move(value);
    return *this;
}

Tag& Tag::operator=(std::vector<int32_t> value)
{
    data_ = std::move(value);
    return *this;
}

Tag& Tag::operator=(std::vector<uint32_t> value)
{
    data_ = std::move(value);
    return *this;
}

Tag& Tag::operator=(std::vector<float> value)
{
    data_ = std::move(value);
    return *this;
}

bool Tag::operator==(const Tag& other) const
{
    return boost::apply_visitor(IsEqualVisitor(), data_, other.data_) &&
           (modifier_ == other.modifier_);
}

bool Tag::operator!=(const Tag& other) const { return !(*this == other); }

bool Tag::HasModifier(const TagModifier m) const
{
    // we just allow one at a time (for now at least)
    return modifier_ == m;
}

bool Tag::IsNull() const { return Type() == TagDataType::INVALID; }

bool Tag::IsInt8() const { return Type() == TagDataType::INT8; }

bool Tag::IsUInt8() const { return Type() == TagDataType::UINT8; }

bool Tag::IsInt16() const { return Type() == TagDataType::INT16; }

bool Tag::IsUInt16() const { return Type() == TagDataType::UINT16; }

bool Tag::IsInt32() const { return Type() == TagDataType::INT32; }

bool Tag::IsUInt32() const { return Type() == TagDataType::UINT32; }

bool Tag::IsFloat() const { return Type() == TagDataType::FLOAT; }

bool Tag::IsString() const { return Type() == TagDataType::STRING; }

bool Tag::IsHexString() const { return IsString() && modifier_ == TagModifier::HEX_STRING; }

bool Tag::IsInt8Array() const { return Type() == TagDataType::INT8_ARRAY; }

bool Tag::IsUInt8Array() const { return Type() == TagDataType::UINT8_ARRAY; }

bool Tag::IsInt16Array() const { return Type() == TagDataType::INT16_ARRAY; }

bool Tag::IsUInt16Array() const { return Type() == TagDataType::UINT16_ARRAY; }

bool Tag::IsInt32Array() const { return Type() == TagDataType::INT32_ARRAY; }

bool Tag::IsUInt32Array() const { return Type() == TagDataType::UINT32_ARRAY; }

bool Tag::IsFloatArray() const { return Type() == TagDataType::FLOAT_ARRAY; }

bool Tag::IsSignedInt() const { return IsInt8() || IsInt16() || IsInt32(); }

bool Tag::IsUnsignedInt() const { return IsUInt8() || IsUInt16() || IsUInt32(); }

bool Tag::IsIntegral() const { return IsSignedInt() || IsUnsignedInt(); }

bool Tag::IsNumeric() const { return IsIntegral() || IsFloat(); }

bool Tag::IsSignedArray() const { return IsInt8Array() || IsInt16Array() || IsInt32Array(); }

bool Tag::IsUnsignedArray() const { return IsUInt8Array() || IsUInt16Array() || IsUInt32Array(); }

bool Tag::IsIntegralArray() const { return IsSignedArray() || IsUnsignedArray(); }

bool Tag::IsArray() const { return IsIntegralArray() || IsFloatArray(); }

TagModifier Tag::Modifier() const { return modifier_; }

Tag& Tag::Modifier(const TagModifier m)
{
    modifier_ = m;
    return *this;
}

char Tag::ToAscii() const { return boost::apply_visitor(AsciiConvertVisitor(), data_); }

int8_t Tag::ToInt8() const
{
    if (IsInt8()) return boost::get<int8_t>(data_);
    return boost::apply_visitor(ToInt8ConvertVisitor(), data_);
}

uint8_t Tag::ToUInt8() const
{
    if (IsUInt8()) return boost::get<uint8_t>(data_);
    return boost::apply_visitor(ToUInt8ConvertVisitor(), data_);
}

int16_t Tag::ToInt16() const
{
    if (IsInt16()) return boost::get<int16_t>(data_);
    return boost::apply_visitor(ToInt16ConvertVisitor(), data_);
}

uint16_t Tag::ToUInt16() const
{
    if (IsUInt16()) return boost::get<uint16_t>(data_);
    return boost::apply_visitor(ToUInt16ConvertVisitor(), data_);
}

int32_t Tag::ToInt32() const
{
    if (IsInt32()) return boost::get<int32_t>(data_);
    return boost::apply_visitor(ToInt32ConvertVisitor(), data_);
}

uint32_t Tag::ToUInt32() const
{
    if (IsUInt32()) return boost::get<uint32_t>(data_);
    return boost::apply_visitor(ToUInt32ConvertVisitor(), data_);
}

float Tag::ToFloat() const { return boost::get<float>(data_); }

std::string Tag::ToString() const { return boost::get<std::string>(data_); }

std::vector<int8_t> Tag::ToInt8Array() const { return boost::get<std::vector<int8_t> >(data_); }

std::vector<uint8_t> Tag::ToUInt8Array() const { return boost::get<std::vector<uint8_t> >(data_); }

std::vector<int16_t> Tag::ToInt16Array() const { return boost::get<std::vector<int16_t> >(data_); }

std::vector<uint16_t> Tag::ToUInt16Array() const
{
    return boost::get<std::vector<uint16_t> >(data_);
}

std::vector<int32_t> Tag::ToInt32Array() const { return boost::get<std::vector<int32_t> >(data_); }

std::vector<uint32_t> Tag::ToUInt32Array() const
{
    return boost::get<std::vector<uint32_t> >(data_);
}

std::vector<float> Tag::ToFloatArray() const { return boost::get<std::vector<float> >(data_); }

TagDataType Tag::Type() const { return TagDataType(data_.which()); }

std::string Tag::Typename() const { return boost::apply_visitor(TypenameVisitor(), data_); }

}  // namespace BAM
}  // namespace PacBio
