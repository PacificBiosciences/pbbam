// File Description
/// \file BamTagCodec.cpp
/// \brief Implements the BamTagCodec class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamTagCodec.h"

#include <cstddef>
#include <cstdint>
#include <cstring>

#include <htslib/kstring.h>

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T>
inline void appendBamValue(const T& value, kstring_t* str)
{
    kputsn_(reinterpret_cast<const char*>(&value), sizeof(value), str);
}

template <typename T>
inline void appendBamMultiValue(const std::vector<T>& container, kstring_t* str)
{
    const uint32_t n = container.size();
    kputsn_(&n, sizeof(n), str);
    kputsn_(reinterpret_cast<const char*>(&container[0]), n * sizeof(T), str);
}

template <typename T>
inline T readBamValue(const uint8_t* src, size_t& offset)
{
    T value;
    memcpy(&value, &src[offset], sizeof(value));
    offset += sizeof(value);
    return value;
}

template <typename T>
std::vector<T> readBamMultiValue(const uint8_t* src, size_t& offset)
{
    uint32_t numElements;
    memcpy(&numElements, &src[offset], sizeof(uint32_t));
    offset += 4;

    std::vector<T> result;
    result.reserve(numElements);
    for (size_t i = 0; i < numElements; ++i) {
        const T value = readBamValue<T>(src, offset);
        result.push_back(value);
    }
    return result;
}

}  // namespace internal

TagCollection BamTagCodec::Decode(const std::vector<uint8_t>& data)
{
    TagCollection tags;

    // NOTE: not completely safe - no real bounds-checking yet on input data

    const uint8_t* pData = data.data();
    const size_t numBytes = data.size();
    size_t i = 0;
    while (i < numBytes) {

        std::string tagName;
        tagName.reserve(2);
        tagName.append(1, pData[i++]);
        tagName.append(1, pData[i++]);

        using internal::readBamMultiValue;
        using internal::readBamValue;

        const auto tagType = static_cast<char>(pData[i++]);
        switch (tagType) {
            case 'A':
            case 'a': {
                tags[tagName] = readBamValue<uint8_t>(pData, i);
                tags[tagName].Modifier(TagModifier::ASCII_CHAR);
                break;
            }

            case 'c':
                tags[tagName] = readBamValue<int8_t>(pData, i);
                break;
            case 'C':
                tags[tagName] = readBamValue<uint8_t>(pData, i);
                break;
            case 's':
                tags[tagName] = readBamValue<int16_t>(pData, i);
                break;
            case 'S':
                tags[tagName] = readBamValue<uint16_t>(pData, i);
                break;
            case 'i':
                tags[tagName] = readBamValue<int32_t>(pData, i);
                break;
            case 'I':
                tags[tagName] = readBamValue<uint32_t>(pData, i);
                break;
            case 'f':
                tags[tagName] = readBamValue<float>(pData, i);
                break;

            case 'Z':
            case 'H': {
                const size_t dataLength = strlen(reinterpret_cast<const char*>(&pData[i]));
                std::string value(reinterpret_cast<const char*>(&pData[i]), dataLength);
                tags[tagName] = value;
                if (tagType == 'H') tags[tagName].Modifier(TagModifier::HEX_STRING);
                i += dataLength + 1;
                break;
            }

            case 'B': {
                const char subTagType = pData[i++];
                switch (subTagType) {
                    case 'c':
                        tags[tagName] = readBamMultiValue<int8_t>(pData, i);
                        break;
                    case 'C':
                        tags[tagName] = readBamMultiValue<uint8_t>(pData, i);
                        break;
                    case 's':
                        tags[tagName] = readBamMultiValue<int16_t>(pData, i);
                        break;
                    case 'S':
                        tags[tagName] = readBamMultiValue<uint16_t>(pData, i);
                        break;
                    case 'i':
                        tags[tagName] = readBamMultiValue<int32_t>(pData, i);
                        break;
                    case 'I':
                        tags[tagName] = readBamMultiValue<uint32_t>(pData, i);
                        break;
                    case 'f':
                        tags[tagName] = readBamMultiValue<float>(pData, i);
                        break;

                    // unknown subTagType
                    default:
                        throw std::runtime_error{"unsupported array-tag-type encountered: " +
                                                 std::string{1, subTagType}};
                }
                break;
            }

            // unknown tagType
            default:
                throw std::runtime_error{"unsupported tag-type encountered: " +
                                         std::string{1, tagType}};
        }
    }

    return tags;
}

std::vector<uint8_t> BamTagCodec::Encode(const TagCollection& tags)
{
    kstring_t str = {0, 0, nullptr};

    for (const auto& tagIter : tags) {

        const auto& name = tagIter.first;
        if (name.size() != 2) throw std::runtime_error{"malformatted tag name: " + name};

        const auto& tag = tagIter.second;
        if (tag.IsNull()) continue;

        // "<TAG>:"
        kputsn_(name.c_str(), 2, &str);

        // "<TYPE>:<DATA>" for printable, ASCII char
        if (tag.HasModifier(TagModifier::ASCII_CHAR)) {
            const char c = tag.ToAscii();
            if (c != '\0') {
                kputc_('A', &str);
                kputc_(c, &str);
                continue;
            }
        }

        using internal::appendBamMultiValue;
        using internal::appendBamValue;

        // "<TYPE>:<DATA>" for all other data
        switch (tag.Type()) {
            case TagDataType::INT8: {
                kputc_('c', &str);
                appendBamValue(tag.ToInt8(), &str);
                break;
            }
            case TagDataType::UINT8: {
                kputc_('C', &str);
                appendBamValue(tag.ToUInt8(), &str);
                break;
            }
            case TagDataType::INT16: {
                kputc_('s', &str);
                appendBamValue(tag.ToInt16(), &str);
                break;
            }
            case TagDataType::UINT16: {
                kputc_('S', &str);
                appendBamValue(tag.ToUInt16(), &str);
                break;
            }
            case TagDataType::INT32: {
                kputc_('i', &str);
                appendBamValue(tag.ToInt32(), &str);
                break;
            }
            case TagDataType::UINT32: {
                kputc_('I', &str);
                appendBamValue(tag.ToUInt32(), &str);
                break;
            }
            case TagDataType::FLOAT: {
                kputc_('f', &str);
                appendBamValue(tag.ToFloat(), &str);
                break;
            }

            case TagDataType::STRING: {
                if (tag.HasModifier(TagModifier::HEX_STRING))
                    kputc_('H', &str);
                else
                    kputc_('Z', &str);
                const auto s = tag.ToString();
                kputsn_(s.c_str(), s.size() + 1, &str);  // this adds the null-term
                break;
            }

            case TagDataType::INT8_ARRAY: {
                kputc_('B', &str);
                kputc_('c', &str);
                appendBamMultiValue(tag.ToInt8Array(), &str);
                break;
            }
            case TagDataType::UINT8_ARRAY: {
                kputc_('B', &str);
                kputc_('C', &str);
                appendBamMultiValue(tag.ToUInt8Array(), &str);
                break;
            }
            case TagDataType::INT16_ARRAY: {
                kputc_('B', &str);
                kputc_('s', &str);
                appendBamMultiValue(tag.ToInt16Array(), &str);
                break;
            }
            case TagDataType::UINT16_ARRAY: {
                kputc_('B', &str);
                kputc_('S', &str);
                appendBamMultiValue(tag.ToUInt16Array(), &str);
                break;
            }
            case TagDataType::INT32_ARRAY: {
                kputc_('B', &str);
                kputc_('i', &str);
                appendBamMultiValue(tag.ToInt32Array(), &str);
                break;
            }
            case TagDataType::UINT32_ARRAY: {
                kputc_('B', &str);
                kputc_('I', &str);
                appendBamMultiValue(tag.ToUInt32Array(), &str);
                break;
            }
            case TagDataType::FLOAT_ARRAY: {
                kputc_('B', &str);
                kputc_('f', &str);
                appendBamMultiValue(tag.ToFloatArray(), &str);
                break;
            }

            // unsupported tag type
            default: {
                free(str.s);
                throw std::runtime_error{"unsupported tag-type encountered: " +
                                         std::to_string(static_cast<uint16_t>(tag.Type()))};
            }
        }
    }

    std::vector<uint8_t> result;
    result.resize(str.l);
    memcpy(reinterpret_cast<char*>(result.data()), str.s, str.l);
    free(str.s);
    return result;
}

Tag BamTagCodec::FromRawData(uint8_t* rawData)
{
    using internal::readBamMultiValue;
    using internal::readBamValue;

    size_t offset = 0;
    const auto tagType = static_cast<char>(*rawData++);
    switch (tagType) {
        case 'A':
        case 'a': {
            Tag t{readBamValue<uint8_t>(rawData, offset)};
            t.Modifier(TagModifier::ASCII_CHAR);
            return t;
        }

        case 'c':
            return {readBamValue<int8_t>(rawData, offset)};
        case 'C':
            return {readBamValue<uint8_t>(rawData, offset)};
        case 's':
            return {readBamValue<int16_t>(rawData, offset)};
        case 'S':
            return {readBamValue<uint16_t>(rawData, offset)};
        case 'i':
            return {readBamValue<int32_t>(rawData, offset)};
        case 'I':
            return {readBamValue<uint32_t>(rawData, offset)};
        case 'f':
            return {readBamValue<float>(rawData, offset)};

        case 'Z':
        case 'H': {
            const size_t dataLength = strlen(reinterpret_cast<const char*>(&rawData[0]));
            std::string value(reinterpret_cast<const char*>(&rawData[0]), dataLength);
            Tag t{value};
            if (tagType == 'H') t.Modifier(TagModifier::HEX_STRING);
            return t;
        }

        case 'B': {
            const char subTagType = *rawData++;
            switch (subTagType) {

                case 'c':
                    return {readBamMultiValue<int8_t>(rawData, offset)};
                case 'C':
                    return {readBamMultiValue<uint8_t>(rawData, offset)};
                case 's':
                    return {readBamMultiValue<int16_t>(rawData, offset)};
                case 'S':
                    return {readBamMultiValue<uint16_t>(rawData, offset)};
                case 'i':
                    return {readBamMultiValue<int32_t>(rawData, offset)};
                case 'I':
                    return {readBamMultiValue<uint32_t>(rawData, offset)};
                case 'f':
                    return {readBamMultiValue<float>(rawData, offset)};

                // unknown subTagType
                default:
                    throw std::runtime_error{"unsupported array-tag-type encountered: " +
                                             std::string{1, subTagType}};
            }
            break;
        }

        // unknown tagType
        default:
            throw std::runtime_error{"unsupported tag-type encountered: " +
                                     std::string{1, tagType}};
    }
    return Tag();  // to avoid compiler warning
}

std::vector<uint8_t> BamTagCodec::ToRawData(const Tag& tag, const TagModifier& additionalModifier)
{
    // temp raw data destination (for use with htslib methods)
    kstring_t str = {0, 0, nullptr};

    // "<TYPE>:<DATA>" for printable, ASCII char
    if (tag.HasModifier(TagModifier::ASCII_CHAR) || additionalModifier == TagModifier::ASCII_CHAR) {
        const char c = tag.ToAscii();
        if (c != '\0') kputc_(c, &str);
    }

    // for all others
    else {

        using internal::appendBamMultiValue;
        using internal::appendBamValue;

        switch (tag.Type()) {

            // single, numeric values
            case TagDataType::INT8:
                appendBamValue(tag.ToInt8(), &str);
                break;
            case TagDataType::UINT8:
                appendBamValue(tag.ToUInt8(), &str);
                break;
            case TagDataType::INT16:
                appendBamValue(tag.ToInt16(), &str);
                break;
            case TagDataType::UINT16:
                appendBamValue(tag.ToUInt16(), &str);
                break;
            case TagDataType::INT32:
                appendBamValue(tag.ToInt32(), &str);
                break;
            case TagDataType::UINT32:
                appendBamValue(tag.ToUInt32(), &str);
                break;
            case TagDataType::FLOAT:
                appendBamValue(tag.ToFloat(), &str);
                break;

            // string & hex-string values
            case TagDataType::STRING: {
                const auto s = tag.ToString();
                kputsn_(s.c_str(), s.size() + 1, &str);  // this adds the null-term
                break;
            }

            // array-type values
            case TagDataType::INT8_ARRAY: {
                kputc_('c', &str);
                appendBamMultiValue(tag.ToInt8Array(), &str);
                break;
            }
            case TagDataType::UINT8_ARRAY: {
                kputc_('C', &str);
                appendBamMultiValue(tag.ToUInt8Array(), &str);
                break;
            }
            case TagDataType::INT16_ARRAY: {
                kputc_('s', &str);
                appendBamMultiValue(tag.ToInt16Array(), &str);
                break;
            }
            case TagDataType::UINT16_ARRAY: {
                kputc_('S', &str);
                appendBamMultiValue(tag.ToUInt16Array(), &str);
                break;
            }
            case TagDataType::INT32_ARRAY: {
                kputc_('i', &str);
                appendBamMultiValue(tag.ToInt32Array(), &str);
                break;
            }
            case TagDataType::UINT32_ARRAY: {
                kputc_('I', &str);
                appendBamMultiValue(tag.ToUInt32Array(), &str);
                break;
            }
            case TagDataType::FLOAT_ARRAY: {
                kputc_('f', &str);
                appendBamMultiValue(tag.ToFloatArray(), &str);
                break;
            }

            // unsupported tag type
            default: {
                free(str.s);
                throw std::runtime_error{"unsupported tag-type encountered: " +
                                         std::to_string(static_cast<uint16_t>(tag.Type()))};
            }
        }
    }

    // store temp contents in actual destination
    std::vector<uint8_t> result;
    result.resize(str.l);
    memcpy(reinterpret_cast<char*>(&result[0]), str.s, str.l);
    free(str.s);
    return result;
}

uint8_t BamTagCodec::TagTypeCode(const Tag& tag, const TagModifier& additionalModifier)
{
    if (tag.HasModifier(TagModifier::ASCII_CHAR) || additionalModifier == TagModifier::ASCII_CHAR) {
        int64_t value = 0;
        switch (tag.Type()) {
            case TagDataType::INT8:
                value = static_cast<int64_t>(tag.ToInt8());
                break;
            case TagDataType::UINT8:
                value = static_cast<int64_t>(tag.ToUInt8());
                break;
            case TagDataType::INT16:
                value = static_cast<int64_t>(tag.ToInt16());
                break;
            case TagDataType::UINT16:
                value = static_cast<int64_t>(tag.ToUInt16());
                break;
            case TagDataType::INT32:
                value = static_cast<int64_t>(tag.ToInt32());
                break;
            case TagDataType::UINT32:
                value = static_cast<int64_t>(tag.ToUInt32());
                break;
            default:
                // non integers not allowed
                throw std::runtime_error{"tag-type not convertible to ASCII, tag-type: " +
                                         std::to_string(static_cast<uint16_t>(tag.Type()))};
        }

        // ensure value is in valid ASCII char range
        if (value < 33 || value > 126)
            throw std::runtime_error{"invalid integer value for ASCII char, value: " +
                                     std::to_string(value)};

        return static_cast<uint8_t>('A');
    }

    switch (tag.Type()) {
        case TagDataType::INT8:
            return static_cast<uint8_t>('c');
        case TagDataType::UINT8:
            return static_cast<uint8_t>('C');
        case TagDataType::INT16:
            return static_cast<uint8_t>('s');
        case TagDataType::UINT16:
            return static_cast<uint8_t>('S');
        case TagDataType::INT32:
            return static_cast<uint8_t>('i');
        case TagDataType::UINT32:
            return static_cast<uint8_t>('I');
        case TagDataType::FLOAT:
            return static_cast<uint8_t>('f');

        case TagDataType::STRING: {
            if (tag.HasModifier(TagModifier::HEX_STRING) ||
                additionalModifier == TagModifier::HEX_STRING)
                return static_cast<uint8_t>('H');
            return static_cast<uint8_t>('Z');
        }

        case TagDataType::INT8_ARRAY:    // fall through
        case TagDataType::UINT8_ARRAY:   // .
        case TagDataType::INT16_ARRAY:   // .
        case TagDataType::UINT16_ARRAY:  // .
        case TagDataType::INT32_ARRAY:   // .
        case TagDataType::UINT32_ARRAY:  // .
        case TagDataType::FLOAT_ARRAY:
            return static_cast<uint8_t>('B');

        default:
            throw std::runtime_error{"unsupported tag-type encountered: " +
                                     std::to_string(static_cast<uint16_t>(tag.Type()))};
    }
    return 0;  // to avoid compiler warning
}

}  // namespace BAM
}  // namespace PacBio
