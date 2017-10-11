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
/// \file SamTagCodec.h
/// \brief Implements the SamTagCodec class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include <cstdint>
#include <limits>

#include "pbbam/SamTagCodec.h"
#include <boost/lexical_cast.hpp>

namespace PacBio {
namespace BAM {
namespace internal {

template<typename T>
inline void appendSamValue(const T& value,
                           std::string& result,
                           bool force8BitInt = false)
{
    if (force8BitInt)
        result.append(boost::lexical_cast<std::string>(static_cast<int>(value)));
    else
        result.append(boost::lexical_cast<std::string>(value));
}

template<typename T>
void appendSamMultiValue(const T& container,
                         std::string& result,
                         bool force8BitInt = false)
{
    auto end = container.cend();
    for (auto iter = container.cbegin(); iter != end; ++iter) {
        result.append(1, ',');
        if ( force8BitInt )
            result.append(boost::lexical_cast<std::string>(static_cast<int>(*iter)));
        else
            result.append(boost::lexical_cast<std::string>(*iter));
    }
}

static
std::vector<std::string> split(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    std::stringstream ss(s);
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
        result.emplace_back(strtof(c+1, &c));
    return result;
}

template<typename T>
std::vector<T> readSignedSamMultiValue(const std::string& data)
{
    std::vector<T> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end)
        result.emplace_back(strtol(c + 1, &c, 0));
    return result;
}

template<typename T>
std::vector<T> readUnsignedSamMultiValue(const std::string& data)
{
    std::vector<T> result;
    auto* c = const_cast<char*>(data.c_str());
    const char* end = c + data.length();
    while (c + 1 < end)
        result.emplace_back(strtoul(c + 1, &c, 0));
    return result;
}

} // namespace internal

TagCollection SamTagCodec::Decode(const std::string& tagString)
{
    TagCollection tags;

    const auto tokens = internal::split(tagString, '\t');
    for (const auto& token : tokens) {
        if (token.size() < 6)          // TT:t:X
            continue;

        const auto name = token.substr(0, 2);
        const auto type = token.at(3);
        const auto remainder = token.substr(5);
        if (remainder.empty())
            throw std::runtime_error("malformatted tag: " + token);

        switch (type) {

            // technically only 'A' is allowed in SAM chars,
            // but we'll be a little permissive
            case 'A' :
            case 'a' :
            {
                tags[name] = Tag(static_cast<char>(remainder.at(0), TagModifier::ASCII_CHAR));
                break;
            }

            // technically only 'i' is allowed in SAM ints, but we'll be a little
            // permissive since SAM might be a bit more "user-edited" than BAM
            case 'c' :
            case 'C' :
            case 's' :
            case 'S' :
            case 'i' :
            case 'I' :
            {
                // check out boost::numeric cast for these conversions

                // negative value (force signed int)
                if (remainder.at(0) == '-') {
                    const auto x = boost::lexical_cast<int32_t>(remainder);
                    if ( x >= std::numeric_limits<int8_t>::min() )
                        tags[name] = static_cast<int8_t>(x);
                    else if ( x >= std::numeric_limits<int16_t>::min() )
                        tags[name] = static_cast<int16_t>(x);
                    else
                        tags[name] = x;
                }

                // unsigned int
                else {
                    const auto x = boost::lexical_cast<uint32_t>(remainder);
                    if ( x <= std::numeric_limits<uint8_t>::max() )
                        tags[name] = static_cast<uint8_t>(x);
                    else if ( x <= std::numeric_limits<uint16_t>::max() )
                        tags[name] = static_cast<uint16_t>(x);
                    else
                        tags[name] = x;
                }
                break;
            }

            case 'f' :
            {
                tags[name] = boost::lexical_cast<float>(remainder);
                break;
            }

            case 'Z' :
            {
                tags[name] = remainder;
                break;
            }

            case 'H' :
            {
                tags[name] = Tag(remainder, TagModifier::HEX_STRING);
                break;
            }

            case 'B' :
            {
                const auto elementType = remainder.at(0);
                const auto arrayData = remainder.substr(1);
                switch (elementType) {
                    case 'c' : tags[name] = internal::readSignedSamMultiValue<int8_t>(arrayData);     break;
                    case 'C' : tags[name] = internal::readUnsignedSamMultiValue<uint8_t>(arrayData);  break;
                    case 's' : tags[name] = internal::readSignedSamMultiValue<int16_t>(arrayData);    break;
                    case 'S' : tags[name] = internal::readUnsignedSamMultiValue<uint16_t>(arrayData); break;
                    case 'i' : tags[name] = internal::readSignedSamMultiValue<int32_t>(arrayData);    break;
                    case 'I' : tags[name] = internal::readUnsignedSamMultiValue<uint32_t>(arrayData); break;
                    case 'f' : tags[name] = internal::readFloatSamMultiValue(arrayData);              break;
                    default:
                        throw std::runtime_error("unsupported array-tag-type encountered: " + std::string(1, elementType));
                }
                break;
            }

            // unsupported SAM tag type
            default :
                throw std::runtime_error("unsupported tag-type encountered: " + std::string(1, type));
        }
    }

    return tags;
}

std::string SamTagCodec::Encode(const TagCollection& tags)
{
    std::string result;
    result.reserve(1024);

    for (const auto& tagIter : tags)
    {
        const auto& name = tagIter.first;
        if (name.size() != 2)
            throw std::runtime_error("malformatted tag name: " + name);

        const auto& tag = tagIter.second;
        if (tag.IsNull())
            continue;

        // tab separator
        if (!result.empty())
            result.append(1, '\t');

        // "<TAG>:"
        result.append(name);
        result.append(1, ':');

        // "<TYPE>:<DATA>" for printable, ASCII char
        if (tag.HasModifier(TagModifier::ASCII_CHAR)) {
            const auto c = tag.ToAscii();
            if (c != '\0') {
                result.append("A:");
                result.append(1, c);
                continue;
            }
        }

        // "<TYPE>:<DATA>" for all other data

        using internal::appendSamMultiValue;
        using internal::appendSamValue;

        switch (tag.Type()) {
            case TagDataType::INT8   : result.append("i:"); appendSamValue(tag.ToInt8(),   result, true); break;
            case TagDataType::UINT8  : result.append("i:"); appendSamValue(tag.ToUInt8(),  result, true); break;
            case TagDataType::INT16  : result.append("i:"); appendSamValue(tag.ToInt16(),  result); break;
            case TagDataType::UINT16 : result.append("i:"); appendSamValue(tag.ToUInt16(), result); break;
            case TagDataType::INT32  : result.append("i:"); appendSamValue(tag.ToInt32(),  result); break;
            case TagDataType::UINT32 : result.append("i:"); appendSamValue(tag.ToUInt32(), result); break;
            case TagDataType::FLOAT  : result.append("f:"); appendSamValue(tag.ToFloat(),  result); break;

            case TagDataType::STRING :
            {
                result.append(tag.HasModifier(TagModifier::HEX_STRING) ? "H:" : "Z:");
                result.append(tag.ToString());
                break;
            }

            case TagDataType::INT8_ARRAY   : result.append("B:c"); appendSamMultiValue(tag.ToInt8Array(),   result, true); break;
            case TagDataType::UINT8_ARRAY  : result.append("B:C"); appendSamMultiValue(tag.ToUInt8Array(),  result, true); break;
            case TagDataType::INT16_ARRAY  : result.append("B:s"); appendSamMultiValue(tag.ToInt16Array(),  result); break;
            case TagDataType::UINT16_ARRAY : result.append("B:S"); appendSamMultiValue(tag.ToUInt16Array(), result); break;
            case TagDataType::INT32_ARRAY  : result.append("B:i"); appendSamMultiValue(tag.ToInt32Array(),  result); break;
            case TagDataType::UINT32_ARRAY : result.append("B:I"); appendSamMultiValue(tag.ToUInt32Array(), result); break;
            case TagDataType::FLOAT_ARRAY  : result.append("B:f"); appendSamMultiValue(tag.ToFloatArray(),  result); break;

            default :
                throw std::runtime_error("unsupported tag-type encountered: " +
                                         std::to_string(static_cast<uint16_t>(tag.Type())));
        }
    }

    return result;
}

} // namespace BAM
} // namespace PacBio
