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
/// \file BamTagCodec.cpp
/// \brief Implements the BamTagCodec class.
//
// Author: Derek Barnett

#include "pbbam/BamTagCodec.h"
#include "AssertUtils.h"
#include <htslib/kstring.h>
#include <cstring>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

template<typename T>
inline void appendBamValue(const T& value, kstring_t* str)
{
    kputsn_((char*)&value, sizeof(value), str);
}

template<typename T>
inline void appendBamMultiValue(const vector<T>& container, kstring_t* str)
{
    const uint32_t n = container.size();
    kputsn_(&n, sizeof(n), str);
    kputsn_((char*)&container[0], n*sizeof(T), str);
}

template<typename T>
inline T readBamValue(const uint8_t* src, size_t& offset)
{
    T value;
    memcpy(&value, &src[offset], sizeof(value));
    offset += sizeof(value);
    return value;
}

template<typename T>
vector<T> readBamMultiValue(const uint8_t* src, size_t& offset)
{
    uint32_t numElements;
    memcpy(&numElements, &src[offset], sizeof(uint32_t));
    offset += 4;

    vector<T> result;
    result.reserve(numElements);
    for (size_t i = 0; i < numElements; ++i) {
        const T& value = readBamValue<T>(src, offset);
        result.push_back(value);
    }
    return result;
}

TagCollection BamTagCodec::Decode(const vector<uint8_t>& data)
{
    TagCollection tags;

    // NOTE: not completely safe - no real bounds-checking yet on input data

    const uint8_t* pData  = data.data();
    const size_t numBytes = data.size();
    size_t i = 0;
    while (i < numBytes) {

        string tagName;
        tagName.reserve(2);
        tagName.append(1, pData[i++]);
        tagName.append(1, pData[i++]);

        const char tagType = static_cast<char>(pData[i++]);
        switch (tagType) {
            case 'A' :
            case 'a' :
            {
                tags[tagName] = readBamValue<uint8_t>(pData, i);
                tags[tagName].Modifier(TagModifier::ASCII_CHAR);
                break;
            }

            case 'c' : tags[tagName] = readBamValue<int8_t>(pData, i);   break;
            case 'C' : tags[tagName] = readBamValue<uint8_t>(pData, i);  break;
            case 's' : tags[tagName] = readBamValue<int16_t>(pData, i);  break;
            case 'S' : tags[tagName] = readBamValue<uint16_t>(pData, i); break;
            case 'i' : tags[tagName] = readBamValue<int32_t>(pData, i);  break;
            case 'I' : tags[tagName] = readBamValue<uint32_t>(pData, i); break;
            case 'f' : tags[tagName] = readBamValue<float>(pData, i);    break;

            case 'Z' :
            case 'H' :
            {
                const size_t dataLength = strlen((const char*)&pData[i]);
                string value;
                value.resize(dataLength);
                memcpy((char*)value.data(), &pData[i], dataLength);
                tags[tagName] = value;
                if (tagType == 'H')
                    tags[tagName].Modifier(TagModifier::HEX_STRING);
                i += dataLength + 1;
                break;
            }

            case 'B' :
            {
                const char subTagType = pData[i++];
                switch (subTagType) {
                    case 'c' : tags[tagName] = readBamMultiValue<int8_t>(pData, i);   break;
                    case 'C' : tags[tagName] = readBamMultiValue<uint8_t>(pData, i);  break;
                    case 's' : tags[tagName] = readBamMultiValue<int16_t>(pData, i);  break;
                    case 'S' : tags[tagName] = readBamMultiValue<uint16_t>(pData, i); break;
                    case 'i' : tags[tagName] = readBamMultiValue<int32_t>(pData, i);  break;
                    case 'I' : tags[tagName] = readBamMultiValue<uint32_t>(pData, i); break;
                    case 'f' : tags[tagName] = readBamMultiValue<float>(pData, i);    break;

                    // unknown subTagType
                    default:
                        PB_ASSERT_OR_RETURN_VALUE(false, TagCollection());
                }
                break;
            }

            // unknown tagType
            default:
                PB_ASSERT_OR_RETURN_VALUE(false, TagCollection());
        }
    }

    return tags;
}

vector<uint8_t> BamTagCodec::Encode(const TagCollection& tags)
{
    kstring_t str = { 0, 0, NULL };

    const auto tagEnd  = tags.cend();
    for (auto tagIter = tags.cbegin(); tagIter != tagEnd; ++tagIter) {
        const string& name = (*tagIter).first;
        const Tag& tag = (*tagIter).second;
        PB_ASSERT_OR_CONTINUE(name.size() == 2);
        if (tag.IsNull())
            continue;

        // "<TAG>:"
        kputsn_(name.c_str(), 2, &str);

        // "<TYPE>:<DATA>" for printable, ASCII char
        if (tag.HasModifier(TagModifier::ASCII_CHAR)) {
            char c = tag.ToAscii();
            if (c != '\0') {
                kputc_('A', &str);
                kputc_(c, &str);
                continue;
            }
        }

        // "<TYPE>:<DATA>" for all other data
        switch ( tag.Type() ) {
            case TagDataType::INT8   :
            {
                kputc_('c', &str);
                appendBamValue(tag.ToInt8(), &str);
                break;
            }
            case TagDataType::UINT8  :
            {
                kputc_('C', &str);
                appendBamValue(tag.ToUInt8(), &str);
                break;
            }
            case TagDataType::INT16  :
            {
                kputc_('s', &str);
                appendBamValue(tag.ToInt16(), &str);
                break;
            }
            case TagDataType::UINT16 :
            {
                kputc_('S', &str);
                appendBamValue(tag.ToUInt16(), &str);
                break;
            }
            case TagDataType::INT32  :
            {
                kputc_('i', &str);
                appendBamValue(tag.ToInt32(), &str);
                break;
            }
            case TagDataType::UINT32 :
            {
                kputc_('I', &str);
                appendBamValue(tag.ToUInt32(), &str);
                break;
            }
            case TagDataType::FLOAT :
            {
                kputc_('f', &str);
                appendBamValue(tag.ToFloat(), &str);
                break;
            }

            case TagDataType::STRING :
            {
                if (tag.HasModifier(TagModifier::HEX_STRING))
                    kputc_('H', &str);
                else
                    kputc_('Z', &str);
                const string& s = tag.ToString();
                kputsn_(s.c_str(), s.size()+1, &str); // this adds the null-term
                break;
            }

            case TagDataType::INT8_ARRAY   :
            {
                kputc_('B', &str);
                kputc_('c', &str);
                appendBamMultiValue(tag.ToInt8Array(), &str);
                break;
            }
            case TagDataType::UINT8_ARRAY  :
            {
                kputc_('B', &str);
                kputc_('C', &str);
                appendBamMultiValue(tag.ToUInt8Array(), &str);
                break;
            }
            case TagDataType::INT16_ARRAY  :
            {
                kputc_('B', &str);
                kputc_('s', &str);
                appendBamMultiValue(tag.ToInt16Array(), &str);
                break;
            }
            case TagDataType::UINT16_ARRAY :
            {
                kputc_('B', &str);
                kputc_('S', &str);
                appendBamMultiValue(tag.ToUInt16Array(), &str);
                break;
            }
            case TagDataType::INT32_ARRAY  :
            {
                kputc_('B', &str);
                kputc_('i', &str);
                appendBamMultiValue(tag.ToInt32Array(), &str);
                break;
            }
            case TagDataType::UINT32_ARRAY :
            {
                kputc_('B', &str);
                kputc_('I', &str);
                appendBamMultiValue(tag.ToUInt32Array(), &str);
                break;
            }
            case TagDataType::FLOAT_ARRAY :
            {
                kputc_('B', &str);
                kputc_('f', &str);
                appendBamMultiValue(tag.ToFloatArray(), &str);
                break;
            }

            // unsupported tag type
            default :
                free(str.s);
                PB_ASSERT_OR_RETURN_VALUE(false, vector<uint8_t>());
        }
    }

    vector<uint8_t> result;
    result.resize(str.l);
    memcpy((char*)&result[0], str.s, str.l);
    free(str.s);
    return result;
}

Tag BamTagCodec::FromRawData(uint8_t* rawData)
{
    size_t offset = 0;
    const char tagType = static_cast<char>(*rawData++);
    switch (tagType) {
        case 'A' :
        case 'a' :
        {
            Tag t = Tag(readBamValue<uint8_t>(rawData, offset));
            t.Modifier(TagModifier::ASCII_CHAR);
            return t;
        }

        case 'c' : return Tag(readBamValue<int8_t>(rawData, offset));
        case 'C' : return Tag(readBamValue<uint8_t>(rawData, offset));
        case 's' : return Tag(readBamValue<int16_t>(rawData, offset));
        case 'S' : return Tag(readBamValue<uint16_t>(rawData, offset));
        case 'i' : return Tag(readBamValue<int32_t>(rawData, offset));
        case 'I' : return Tag(readBamValue<uint32_t>(rawData, offset));
        case 'f' : return Tag(readBamValue<float>(rawData, offset));

        case 'Z' :
        case 'H' :
        {
            const size_t dataLength = strlen((const char*)&rawData[0]);
            string value;
            value.resize(dataLength);
            memcpy((char*)value.data(), &rawData[0], dataLength);
            Tag t(value);
            if (tagType == 'H')
                t.Modifier(TagModifier::HEX_STRING);
            return t;
        }

        case 'B' :
        {
            const char subTagType = *rawData++;
            switch (subTagType) {

                case 'c' : return Tag(readBamMultiValue<int8_t>(rawData, offset));
                case 'C' : return Tag(readBamMultiValue<uint8_t>(rawData, offset));
                case 's' : return Tag(readBamMultiValue<int16_t>(rawData, offset));
                case 'S' : return Tag(readBamMultiValue<uint16_t>(rawData, offset));
                case 'i' : return Tag(readBamMultiValue<int32_t>(rawData, offset));
                case 'I' : return Tag(readBamMultiValue<uint32_t>(rawData, offset));
                case 'f' : return Tag(readBamMultiValue<float>(rawData, offset));

                // unknown subTagType
                default:
                    PB_ASSERT_OR_RETURN_VALUE(false, Tag());
            }
            break;
        }

        // unknown tagType
        default:
            PB_ASSERT_OR_RETURN_VALUE(false, Tag());
    }
    return Tag(); // to avoid compiler warning
}

vector<uint8_t> BamTagCodec::ToRawData(const Tag& tag,
                                       const TagModifier& additionalModifier)
{
    // temp raw data destination (for use with htslib methods)
    kstring_t str = { 0, 0, NULL };

    // "<TYPE>:<DATA>" for printable, ASCII char
    if (tag.HasModifier(TagModifier::ASCII_CHAR) || additionalModifier == TagModifier::ASCII_CHAR) {
        const char c = tag.ToAscii();
        if (c != '\0')
            kputc_(c, &str);
    }

    // for all others
    else {
        switch (tag.Type()) {

            // single, numeric values
            case TagDataType::INT8   : appendBamValue(tag.ToInt8(),   &str); break;
            case TagDataType::UINT8  : appendBamValue(tag.ToUInt8(),  &str); break;
            case TagDataType::INT16  : appendBamValue(tag.ToInt16(),  &str); break;
            case TagDataType::UINT16 : appendBamValue(tag.ToUInt16(), &str); break;
            case TagDataType::INT32  : appendBamValue(tag.ToInt32(),  &str); break;
            case TagDataType::UINT32 : appendBamValue(tag.ToUInt32(), &str); break;
            case TagDataType::FLOAT  : appendBamValue(tag.ToFloat(),  &str); break;

            // string & hex-string values
            case TagDataType::STRING :
            {
                const string& s = tag.ToString();
                kputsn_(s.c_str(), s.size()+1, &str); // this adds the null-term
                break;
            }

            // array-type values
            case TagDataType::INT8_ARRAY   :
            {
                kputc_('c', &str);
                appendBamMultiValue(tag.ToInt8Array(), &str);
                break;
            }
            case TagDataType::UINT8_ARRAY  :
            {
                kputc_('C', &str);
                appendBamMultiValue(tag.ToUInt8Array(), &str);
                break;
            }
            case TagDataType::INT16_ARRAY  :
            {
                kputc_('s', &str);
                appendBamMultiValue(tag.ToInt16Array(), &str);
                break;
            }
            case TagDataType::UINT16_ARRAY :
            {
                kputc_('S', &str);
                appendBamMultiValue(tag.ToUInt16Array(), &str);
                break;
            }
            case TagDataType::INT32_ARRAY  :
            {
                kputc_('i', &str);
                appendBamMultiValue(tag.ToInt32Array(), &str);
                break;
            }
            case TagDataType::UINT32_ARRAY :
            {
                kputc_('I', &str);
                appendBamMultiValue(tag.ToUInt32Array(), &str);
                break;
            }
            case TagDataType::FLOAT_ARRAY :
            {
                kputc_('f', &str);
                appendBamMultiValue(tag.ToFloatArray(), &str);
                break;
            }

            // unsupported tag type
            default :
                free(str.s);
                PB_ASSERT_OR_RETURN_VALUE(false, vector<uint8_t>());
        }
    }

    // store temp contents in actual destination
    vector<uint8_t> result;
    result.resize(str.l);
    memcpy((char*)&result[0], str.s, str.l);
    free(str.s);
    return result;
}

uint8_t BamTagCodec::TagTypeCode(const Tag& tag,
                                 const TagModifier& additionalModifier)
{
    if (tag.HasModifier(TagModifier::ASCII_CHAR) || additionalModifier == TagModifier::ASCII_CHAR) {
        int64_t value = 0;
        switch (tag.Type()) {
            case TagDataType::INT8   : value = static_cast<int64_t>(tag.ToInt8());   break;
            case TagDataType::UINT8  : value = static_cast<int64_t>(tag.ToUInt8());  break;
            case TagDataType::INT16  : value = static_cast<int64_t>(tag.ToInt16());  break;
            case TagDataType::UINT16 : value = static_cast<int64_t>(tag.ToUInt16()); break;
            case TagDataType::INT32  : value = static_cast<int64_t>(tag.ToInt32());  break;
            case TagDataType::UINT32 : value = static_cast<int64_t>(tag.ToUInt32()); break;
            default:
                // non integers not allowed
                PB_ASSERT_OR_RETURN_VALUE(false, 0);
        }
        // printable range
        PB_ASSERT_OR_RETURN_VALUE(value >= 33,  0);
        PB_ASSERT_OR_RETURN_VALUE(value <= 126, 0);
        return static_cast<uint8_t>('A');
    }

    switch (tag.Type()) {
        case TagDataType::INT8   : return static_cast<uint8_t>('c');
        case TagDataType::UINT8  : return static_cast<uint8_t>('C');
        case TagDataType::INT16  : return static_cast<uint8_t>('s');
        case TagDataType::UINT16 : return static_cast<uint8_t>('S');
        case TagDataType::INT32  : return static_cast<uint8_t>('i');
        case TagDataType::UINT32 : return static_cast<uint8_t>('I');
        case TagDataType::FLOAT  : return static_cast<uint8_t>('f');

        case TagDataType::STRING :
        {
            if (tag.HasModifier(TagModifier::HEX_STRING) || additionalModifier == TagModifier::HEX_STRING)
                return static_cast<uint8_t>('H');
            else
                return static_cast<uint8_t>('Z');
        }

        case TagDataType::INT8_ARRAY   : // fall through
        case TagDataType::UINT8_ARRAY  : // .
        case TagDataType::INT16_ARRAY  : // .
        case TagDataType::UINT16_ARRAY : // .
        case TagDataType::INT32_ARRAY  : // .
        case TagDataType::UINT32_ARRAY : // .
        case TagDataType::FLOAT_ARRAY  : return static_cast<uint8_t>('B');

        default:
            PB_ASSERT_OR_RETURN_VALUE(false, 0);
    }
    return 0; // to avoid compiler warning
}
