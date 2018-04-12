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

Tag::Tag(int8_t value) : data_(value) {}
Tag::Tag(uint8_t value) : data_(value) {}
Tag::Tag(int16_t value) : data_(value) {}
Tag::Tag(uint16_t value) : data_(value) {}
Tag::Tag(int32_t value) : data_(value) {}
Tag::Tag(uint32_t value) : data_(value) {}
Tag::Tag(float value) : data_(value) {}
Tag::Tag(const std::string& value) : data_(value) {}
Tag::Tag(const std::vector<int8_t>& value) : data_(value) {}
Tag::Tag(const std::vector<uint8_t>& value) : data_(value) {}
Tag::Tag(const std::vector<int16_t>& value) : data_(value) {}
Tag::Tag(const std::vector<uint16_t>& value) : data_(value) {}
Tag::Tag(const std::vector<int32_t>& value) : data_(value) {}
Tag::Tag(const std::vector<uint32_t>& value) : data_(value) {}
Tag::Tag(const std::vector<float>& value) : data_(value) {}

Tag::Tag(int8_t value, const TagModifier mod) : data_(value), modifier_(mod)
{
    if (mod == TagModifier::HEX_STRING)
        throw std::runtime_error(
            "HEX_STRING is not a valid tag modifier for int8_t data. "
            "It is intended for string-type data only.");
}

Tag::Tag(const std::string& value, const TagModifier mod) : data_(value), modifier_(mod)
{
    if (mod == TagModifier::ASCII_CHAR)
        throw std::runtime_error(
            "ASCII_CHAR is not a valid tag modifier for string-type data. "
            "To construct an ASCII char tag, use a single-quoted value (e.g. 'X' instead of "
            "\"X\")");
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
Tag& Tag::operator=(const std::string& value)
{
    data_ = value;
    return *this;
}
Tag& Tag::operator=(const std::vector<int8_t>& value)
{
    data_ = value;
    return *this;
}
Tag& Tag::operator=(const std::vector<uint8_t>& value)
{
    data_ = value;
    return *this;
}
Tag& Tag::operator=(const std::vector<int16_t>& value)
{
    data_ = value;
    return *this;
}
Tag& Tag::operator=(const std::vector<uint16_t>& value)
{
    data_ = value;
    return *this;
}
Tag& Tag::operator=(const std::vector<int32_t>& value)
{
    data_ = value;
    return *this;
}
Tag& Tag::operator=(const std::vector<uint32_t>& value)
{
    data_ = value;
    return *this;
}
Tag& Tag::operator=(const std::vector<float>& value)
{
    data_ = value;
    return *this;
}

}  // namespace BAM
}  // namespace PacBio
