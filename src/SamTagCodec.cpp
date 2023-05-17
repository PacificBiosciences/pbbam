#include "PbbamInternalConfig.h"

#include <pbbam/SamTagCodec.h>

#include <pbbam/StringUtilities.h>

#include <boost/lexical_cast.hpp>

#include <limits>
#include <sstream>
#include <string>

#include <cstdint>

namespace PacBio {
namespace BAM {
namespace {

std::vector<float> readFloatSamMultiValue(const std::string& data)
{
    std::vector<float> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end) {
        result.emplace_back(strtof(c + 1, &c));
    }
    return result;
}

template <typename T>
std::vector<T> readSignedSamMultiValue(const std::string& data)
{
    std::vector<T> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end) {
        result.emplace_back(strtol(c + 1, &c, 0));
    }
    return result;
}

template <typename T>
std::vector<T> readUnsignedSamMultiValue(const std::string& data)
{
    std::vector<T> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end) {
        result.emplace_back(strtoul(c + 1, &c, 0));
    }
    return result;
}

}  // namespace

TagCollection SamTagCodec::Decode(const std::string& tagString)
{
    TagCollection tags;

    const auto tokens = Split(tagString, '\t');
    for (const auto& token : tokens) {
        if (token.size() < 6) {  // TT:t:X
            continue;
        }

        const auto name = token.substr(0, 2);
        const auto type = token.at(3);
        const auto remainder = token.substr(5);
        if (remainder.empty()) {
            throw std::runtime_error{"[pbbam] SAM tag ERROR: malformed tag: " + token};
        }

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
                    const auto x = boost::lexical_cast<std::int32_t>(remainder);
                    if (x >= std::numeric_limits<std::int8_t>::min()) {
                        tags[name] = static_cast<std::int8_t>(x);
                    } else if (x >= std::numeric_limits<std::int16_t>::min()) {
                        tags[name] = static_cast<std::int16_t>(x);
                    } else {
                        tags[name] = x;
                    }
                }

                // unsigned int
                else {
                    const auto x = boost::lexical_cast<std::uint32_t>(remainder);
                    if (x <= std::numeric_limits<std::uint8_t>::max()) {
                        tags[name] = static_cast<std::uint8_t>(x);
                    } else if (x <= std::numeric_limits<std::uint16_t>::max()) {
                        tags[name] = static_cast<std::uint16_t>(x);
                    } else {
                        tags[name] = x;
                    }
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
                        tags[name] = readSignedSamMultiValue<std::int8_t>(arrayData);
                        break;
                    case 'C':
                        tags[name] = readUnsignedSamMultiValue<std::uint8_t>(arrayData);
                        break;
                    case 's':
                        tags[name] = readSignedSamMultiValue<std::int16_t>(arrayData);
                        break;
                    case 'S':
                        tags[name] = readUnsignedSamMultiValue<std::uint16_t>(arrayData);
                        break;
                    case 'i':
                        tags[name] = readSignedSamMultiValue<std::int32_t>(arrayData);
                        break;
                    case 'I':
                        tags[name] = readUnsignedSamMultiValue<std::uint32_t>(arrayData);
                        break;
                    case 'f':
                        tags[name] = readFloatSamMultiValue(arrayData);
                        break;
                    default:
                        throw std::runtime_error{
                            "[pbbam] SAM tag ERROR: unsupported array-tag-type encountered: " +
                            std::string{1, elementType}};
                }
                break;
            }

            // unsupported SAM tag type
            default:
                throw std::runtime_error{
                    "[pbbam] SAM tag ERROR: unsupported tag-type encountered: " +
                    std::string{1, type}};
        }
    }

    return tags;
}

std::string SamTagCodec::Encode(const std::string& name, const PacBio::BAM::Tag& tag)
{
    // upfront checks
    if (name.size() != 2) {
        throw std::runtime_error{"[pbbam] SAM tag ERROR: malformed tag name: " + name};
    }
    if (tag.IsNull()) {
        return {};
    }

    // "<TAG>:"
    std::ostringstream result;
    result << name << ':';

    // ASCII char
    if (tag.HasModifier(TagModifier::ASCII_CHAR)) {
        const auto c = tag.ToAscii();
        if (c != '\0') {
            result << "A:" << c;
            return result.str();
        }
    }

    // "<TYPE>:<DATA>" for all other data
    switch (tag.Type()) {
        case TagDataType::INT8:
            result << "i:" << static_cast<std::int32_t>(tag.ToInt8());
            break;
        case TagDataType::UINT8:
            result << "i:" << static_cast<std::int32_t>(tag.ToUInt8());
            break;
        case TagDataType::INT16:
            result << "i:" << tag.ToInt16();
            break;
        case TagDataType::UINT16:
            result << "i:" << tag.ToUInt16();
            break;
        case TagDataType::INT32:
            result << "i:" << tag.ToInt32();
            break;
        case TagDataType::UINT32:
            result << "i:" << tag.ToUInt32();
            break;
        case TagDataType::FLOAT:
            result << "f:" << tag.ToFloat();
            break;
        case TagDataType::STRING:
            result << (tag.HasModifier(TagModifier::HEX_STRING) ? 'H' : 'Z') << ':'
                   << tag.ToString();
            break;
        case TagDataType::INT8_ARRAY:
            result << "B:c";
            for (const std::int8_t x : tag.ToInt8Array()) {
                result << ',' << static_cast<std::int32_t>(x);
            }
            break;
        case TagDataType::UINT8_ARRAY:
            result << "B:C";
            for (const std::uint8_t x : tag.ToUInt8Array()) {
                result << ',' << static_cast<std::uint32_t>(x);
            }
            break;
        case TagDataType::INT16_ARRAY:
            result << "B:s";
            for (const std::int16_t x : tag.ToInt16Array()) {
                result << ',' << x;
            }
            break;
        case TagDataType::UINT16_ARRAY:
            result << "B:S";
            for (const std::uint16_t x : tag.ToUInt16Array()) {
                result << ',' << x;
            }
            break;
        case TagDataType::INT32_ARRAY:
            result << "B:i";
            for (const std::int32_t x : tag.ToInt32Array()) {
                result << ',' << x;
            }
            break;
        case TagDataType::UINT32_ARRAY:
            result << "B:I";
            for (const std::uint32_t x : tag.ToUInt32Array()) {
                result << ',' << x;
            }
            break;
        case TagDataType::FLOAT_ARRAY:
            result << "B:f";
            for (const float x : tag.ToFloatArray()) {
                result << ',' << x;
            }
            break;
        default:
            throw std::runtime_error{"[pbbam] SAM tag ERROR: unsupported tag-type encountered: " +
                                     std::to_string(static_cast<std::uint16_t>(tag.Type()))};
    }
    return result.str();
}

std::string SamTagCodec::Encode(const TagCollection& tags)
{
    std::ostringstream result;
    for (const auto& tagIter : tags) {
        const std::string& name = tagIter.first;
        const Tag& tag = tagIter.second;
        if (!result.str().empty()) {
            result << '\t';
        }
        result << Encode(name, tag);
    }
    return result.str();
}

std::string MakeSamTag(std::string tag, std::string value)
{
    return '\t' + std::move(tag) + ':' + std::move(value);
}

}  // namespace BAM
}  // namespace PacBio
