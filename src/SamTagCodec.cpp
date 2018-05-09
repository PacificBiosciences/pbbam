// File Description
/// \file SamTagCodec.h
/// \brief Implements the SamTagCodec class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SamTagCodec.h"

#include <cstdint>
#include <limits>

#include <boost/lexical_cast.hpp>

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T>
inline void appendSamValue(const T& value, std::string& result)
{
    result.append(boost::lexical_cast<std::string>(value));
}

template <typename T>
inline void appendSamValue_8bit(const T& value, std::string& result)
{
    result.append(boost::lexical_cast<std::string>(static_cast<int>(value)));
}

template <typename T>
void appendSamMultiValue(const T& container, std::string& result)
{
    for (const auto x : container) {
        result.push_back(',');
        appendSamValue(x, result);
    }
}

template <typename T>
void appendSamMultiValue_8bit(const T& container, std::string& result)
{
    for (const auto x : container) {
        result.push_back(',');
        appendSamValue_8bit(x, result);
    }
}

static std::vector<std::string> split(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    std::istringstream ss{s};
    std::string item;
    while (std::getline(ss, item, delim))
        elems.push_back(item);
    return elems;
}

std::vector<float> readFloatSamMultiValue(const std::string& data)
{
    std::vector<float> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end)
        result.emplace_back(strtof(c + 1, &c));
    return result;
}

template <typename T>
std::vector<T> readSignedSamMultiValue(const std::string& data)
{
    std::vector<T> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end)
        result.emplace_back(strtol(c + 1, &c, 0));
    return result;
}

template <typename T>
std::vector<T> readUnsignedSamMultiValue(const std::string& data)
{
    std::vector<T> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end)
        result.emplace_back(strtoul(c + 1, &c, 0));
    return result;
}

}  // namespace internal

TagCollection SamTagCodec::Decode(const std::string& tagString)
{
    TagCollection tags;

    const auto tokens = internal::split(tagString, '\t');
    for (const auto& token : tokens) {
        if (token.size() < 6)  // TT:t:X
            continue;

        const auto name = token.substr(0, 2);
        const auto type = token.at(3);
        const auto remainder = token.substr(5);
        if (remainder.empty()) throw std::runtime_error{"malformatted tag: " + token};

        switch (type) {

            // technically only 'A' is allowed in SAM chars,
            // but we'll be a little permissive
            case 'A':
            case 'a': {
                tags[name] = Tag{static_cast<char>(remainder[0], TagModifier::ASCII_CHAR)};
                break;
            }

            // technically only 'i' is allowed in SAM ints, but we'll be a little
            // permissive since SAM might be a bit more "user-edited" than BAM
            case 'c':
            case 'C':
            case 's':
            case 'S':
            case 'i':
            case 'I': {
                // check out boost::numeric cast for these conversions

                // negative value (force signed int)
                if (remainder[0] == '-') {
                    const auto x = boost::lexical_cast<int32_t>(remainder);
                    if (x >= std::numeric_limits<int8_t>::min())
                        tags[name] = static_cast<int8_t>(x);
                    else if (x >= std::numeric_limits<int16_t>::min())
                        tags[name] = static_cast<int16_t>(x);
                    else
                        tags[name] = x;
                }

                // unsigned int
                else {
                    const auto x = boost::lexical_cast<uint32_t>(remainder);
                    if (x <= std::numeric_limits<uint8_t>::max())
                        tags[name] = static_cast<uint8_t>(x);
                    else if (x <= std::numeric_limits<uint16_t>::max())
                        tags[name] = static_cast<uint16_t>(x);
                    else
                        tags[name] = x;
                }
                break;
            }

            case 'f': {
                tags[name] = boost::lexical_cast<float>(remainder);
                break;
            }

            case 'Z': {
                tags[name] = remainder;
                break;
            }

            case 'H': {
                tags[name] = Tag(remainder, TagModifier::HEX_STRING);
                break;
            }

            case 'B': {
                const auto elementType = remainder[0];
                const auto arrayData = remainder.substr(1);
                switch (elementType) {
                    case 'c':
                        tags[name] = internal::readSignedSamMultiValue<int8_t>(arrayData);
                        break;
                    case 'C':
                        tags[name] = internal::readUnsignedSamMultiValue<uint8_t>(arrayData);
                        break;
                    case 's':
                        tags[name] = internal::readSignedSamMultiValue<int16_t>(arrayData);
                        break;
                    case 'S':
                        tags[name] = internal::readUnsignedSamMultiValue<uint16_t>(arrayData);
                        break;
                    case 'i':
                        tags[name] = internal::readSignedSamMultiValue<int32_t>(arrayData);
                        break;
                    case 'I':
                        tags[name] = internal::readUnsignedSamMultiValue<uint32_t>(arrayData);
                        break;
                    case 'f':
                        tags[name] = internal::readFloatSamMultiValue(arrayData);
                        break;
                    default:
                        throw std::runtime_error{"unsupported array-tag-type encountered: " +
                                                 std::string{1, elementType}};
                }
                break;
            }

            // unsupported SAM tag type
            default:
                throw std::runtime_error{"unsupported tag-type encountered: " +
                                         std::string{1, type}};
        }
    }

    return tags;
}

std::string SamTagCodec::Encode(const TagCollection& tags)
{
    std::string result;
    result.reserve(1024);

    for (const auto& tagIter : tags) {
        const auto& name = tagIter.first;
        if (name.size() != 2) throw std::runtime_error{"malformatted tag name: " + name};

        const auto& tag = tagIter.second;
        if (tag.IsNull()) continue;

        // tab separator
        if (!result.empty()) result.push_back('\t');

        // "<TAG>:"
        result.append(name);
        result.push_back(':');

        // "<TYPE>:<DATA>" for printable, ASCII char
        if (tag.HasModifier(TagModifier::ASCII_CHAR)) {
            const auto c = tag.ToAscii();
            if (c != '\0') {
                result.push_back('A');
                result.push_back(':');
                result.push_back(c);
                continue;
            }
        }

        // "<TYPE>:<DATA>" for all other data

        using internal::appendSamMultiValue;
        using internal::appendSamMultiValue_8bit;
        using internal::appendSamValue;
        using internal::appendSamValue_8bit;

        switch (tag.Type()) {
            case TagDataType::INT8:
                result.push_back('i');
                result.push_back(':');
                appendSamValue_8bit(tag.ToInt8(), result);
                break;
            case TagDataType::UINT8:
                result.push_back('i');
                result.push_back(':');
                appendSamValue_8bit(tag.ToUInt8(), result);
                break;
            case TagDataType::INT16:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToInt16(), result);
                break;
            case TagDataType::UINT16:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToUInt16(), result);
                break;
            case TagDataType::INT32:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToInt32(), result);
                break;
            case TagDataType::UINT32:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToUInt32(), result);
                break;
            case TagDataType::FLOAT:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToFloat(), result);
                break;

            case TagDataType::STRING: {
                result.push_back(tag.HasModifier(TagModifier::HEX_STRING) ? 'H' : 'Z');
                result.push_back(':');
                result.append(tag.ToString());
                break;
            }

            case TagDataType::INT8_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('c');
                appendSamMultiValue_8bit(tag.ToInt8Array(), result);
                break;
            case TagDataType::UINT8_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('C');
                appendSamMultiValue_8bit(tag.ToUInt8Array(), result);
                break;
            case TagDataType::INT16_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('s');
                appendSamMultiValue(tag.ToInt16Array(), result);
                break;
            case TagDataType::UINT16_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('S');
                appendSamMultiValue(tag.ToUInt16Array(), result);
                break;
            case TagDataType::INT32_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('i');
                appendSamMultiValue(tag.ToInt32Array(), result);
                break;
            case TagDataType::UINT32_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('I');
                appendSamMultiValue(tag.ToUInt32Array(), result);
                break;
            case TagDataType::FLOAT_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('f');
                appendSamMultiValue(tag.ToFloatArray(), result);
                break;

            default:
                throw std::runtime_error{"unsupported tag-type encountered: " +
                                         std::to_string(static_cast<uint16_t>(tag.Type()))};
        }
    }

    return result;
}

}  // namespace BAM
}  // namespace PacBio
