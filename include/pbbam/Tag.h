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
/// \file Tag.h
/// \brief Defines the Tag class.
//
// Author: Derek Barnett

#ifndef TAG_H
#define TAG_H

#include "pbbam/Config.h"
#include <boost/variant.hpp>
#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief This enum is used to describe the exact (C++) data type held by a
///        Tag.
///
enum class TagDataType
{
    INVALID      = 0     ///< boost::blank
  , INT8                 ///< int8_t
  , UINT8                ///< uint8_t
  , INT16                ///< int16_t
  , UINT16               ///< uint16_t
  , INT32        = 5     ///< int32_t
  , UINT32               ///< uint32_t
  , FLOAT                ///< float
  , STRING               ///< std::string
  , INT8_ARRAY           ///< std::vector<int8_t>
  , UINT8_ARRAY  = 10    ///< std::vector<uint8_t>
  , INT16_ARRAY          ///< std::vector<int16_t>
  , UINT16_ARRAY         ///< std::vector<uint16_t>
  , INT32_ARRAY          ///< std::vector<int32_t>
  , UINT32_ARRAY         ///< std::vector<uint32_t>
  , FLOAT_ARRAY  = 15    ///< std::vector<float>
};

/// \brief This enum provides additional instructions on interpreting the tag's
///        value.
///
/// Some C++ data types (e.g. std::string) may represent more than one BAM tag
/// type ('H' vs 'Z'). Thus a TagModifier may be used to indicate how to
/// properly distinguish between these shared data types.
///
enum class TagModifier
{
    /// \brief This value indicates that the tag has no modifiers set.
    ///
    NONE = 0,

    /// \brief This modifier marks an integer as ASCII.
    ///
    /// SAM/BAM has the concept of an ASCII character that is distinct from an
    /// 8-bit integer. However, there is no such pure separation in C++ - as
    /// int8_t/uint8_t are likely implemented as typedefs around char/unsigned
    /// char. Thus this modifier can be used to indicate a tag's value should be
    /// interpreted as a printable, ASCII character.
    ///
    ASCII_CHAR,

    /// \brief This modifier marks std::string data as "hex string", rather than
    ///        a regular string.
    ///
    /// SAM/BAM has a distinction between regular strings and "Hex format"
    /// strings. However, they are both manipulated in C++ via std::string. Thus
    /// this modifier can be used to indicate that a tag's string data should be
    /// interpreted as "Hex format" rather than a regular, literal string.
    ///
    HEX_STRING
};

/// \brief The Tag class represents a SAM/BAM record tag value.
///
/// SAM/BAM tags may store values from a variety of types: varying fixed-width
/// integers, strings, arrays of data, etc.
///
/// The Tag class allow tags to be handled in a generic fashion, while
/// maintaining a high level of type-safety. Only those types recognized by the
/// SAM/BAM format are allowed, and extracting the value from a tag is subject
/// to allowed conversion rules, as well.
///
// Inspired by (but greatly simplified & modified from) the boost::variant
// wrapper approach taken by DynamicCpp (https://code.google.com/p/dynamic-cpp)
//
class PBBAM_EXPORT Tag
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty, null tag
    Tag() = default;

    /// \brief Creates a Tag from a signed 8-bit integer or character.
    ///
    /// Without a TagModifier, the resulting Tag will be annotated as containing
    /// an 8-bit integer, whether the input \p value was an integer or a char.
    /// For ASCII tags, use one of these methods:
    /// \include code/Tag_AsciiCtor.txt
    ///
    Tag(int8_t value);

   /// \brief Creates a Tag from a signed 8-bit integer or character,
   ///        applying the provided modifier.
    ///
    /// This method allows direct construction of an ASCII character, rather
    /// than an 8-bit integer (e.g. Tag('A', TagModifier::ASCII_CHAR) ).
    ///
    /// \throws runtime_error if \p modifier is not valid for int8_t data
    ///
    Tag(int8_t value, const TagModifier mod);

    /// \brief Creates a Tag from an unsigned 8-bit integer or character.
    ///
    /// Without a TagModifier, the resulting Tag will be annotated as containing
    /// an 8-bit unsigned integer, whether the input \p value was an integer or
    /// a char. For ASCII tags, use one of these methods:
    /// \include code/Tag_AsciiCtor.txt
    ///
    Tag(uint8_t value);

    /// \brief Creates a Tag from 16-bit integer.
    Tag(int16_t value);

    /// \brief Creates a Tag from 16-bit unsigned integer.
    Tag(uint16_t value);

    /// \brief Creates a Tag from 32-bit signed integer.
    Tag(int32_t value);

    /// \brief Creates a Tag from 32-bit unsigned integer.
    Tag(uint32_t value);

    /// \brief Creates a Tag from floating-point value.
    Tag(float value);

    /// \brief Creates a Tag from string data.
    Tag(const std::string& value);

    /// \brief Creates a Tag from string data, adding modifier.
    ///
    /// \throws runtime_error if \p modifier is not valid for string data
    ///
    Tag(const std::string& value, const TagModifier mod);

    /// \brief Creates a Tag from a vector of 8-bit integers.
    Tag(const std::vector<int8_t>& value);

    /// \brief Creates a Tag from a vector of 8-bit unsigned integers.
    Tag(const std::vector<uint8_t>& value);

    /// \brief Creates a Tag from a vector of 16-bit integers.
    Tag(const std::vector<int16_t>& value);

    /// \brief Creates a Tag from a vector of 16-bit unsigned integers.
    Tag(const std::vector<uint16_t>& value);

    /// Constructs a Tag from a vector of 32-bit integers.
    Tag(const std::vector<int32_t>& value);

    /// \brief Creates a Tag from a vector of 32-bit unsigned integers.
    Tag(const std::vector<uint32_t>& value);

    /// \brief Creates a Tag from a vector of floating-point values.
    Tag(const std::vector<float>& value);
    
    Tag(const Tag& other) = default;
    Tag(Tag&& other) = default;
    ~Tag() = default;

    Tag& operator=(boost::blank value);
    Tag& operator=(int8_t value);
    Tag& operator=(uint8_t value);
    Tag& operator=(int16_t value);
    Tag& operator=(uint16_t value);
    Tag& operator=(int32_t value);
    Tag& operator=(uint32_t value);
    Tag& operator=(float value);
    Tag& operator=(const std::string& value);
    Tag& operator=(const std::vector<int8_t>& value);
    Tag& operator=(const std::vector<uint8_t>& value);
    Tag& operator=(const std::vector<int16_t>& value);
    Tag& operator=(const std::vector<uint16_t>& value);
    Tag& operator=(const std::vector<int32_t>& value);
    Tag& operator=(const std::vector<uint32_t>& value);
    Tag& operator=(const std::vector<float>& value);

    Tag& operator=(const Tag& other) = default;
    Tag& operator=(Tag&& other) = default;

    bool operator== (const Tag& other) const;
    bool operator!= (const Tag& other) const;

    /// \}

public:
    /// \name Data Conversion & Validation
    /// \{

    /// \brief Converts the tag value to an ASCII character.
    ///
    /// Tag must hold an integral type, within the valid ASCII range [33-127].
    ///
    /// \returns ASCII character value
    /// \throws std::runtime_error if not ASCII-compatible
    ///
    char ToAscii() const;

    /// \returns tag data as signed 8-bit (casting if needed)
    /// \throws std::runtime_error if not integral data, or out of valid range
    int8_t ToInt8() const;

    /// \returns tag data as unsigned 8-bit (casting if needed)
    /// \throws std::runtime_error if not integral data, or out of valid range
    uint8_t ToUInt8() const;

    /// \returns tag data as signed 16-bit (casting if needed)
    /// \throws std::runtime_error if not integral data, or out of valid range
    int16_t ToInt16() const;

    /// \returns tag data as unsigned 16-bit (casting if needed)
    /// \throws std::runtime_error if not integral data, or out of valid range
    uint16_t ToUInt16() const;

    /// \returns tag data as signed 32-bit (casting if needed)
    /// \throws std::runtime_error if not integral data, or out of valid range
    int32_t ToInt32() const;

    /// \returns tag data as unsigned 32-bit (casting if needed)
    /// \throws std::runtime_error if not integral data, or out of valid range
    uint32_t ToUInt32() const;

    /// \returns tag data as float
    /// \throws std::runtime_error if tag does not contain a value of
    ///         explicit type: float
    float ToFloat() const;

    /// \returns tag data as std::string
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::string
    std::string ToString() const;

    /// \returns tag data as std::vector<int8_t>
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::vector<int8_t>
    std::vector<int8_t> ToInt8Array() const;

    /// \returns tag data as std::vector<uint8_t>
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::vector<uint8_t>
    std::vector<uint8_t> ToUInt8Array() const;

    /// \returns tag data as std::vector<int16_t>
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::vector<int16_t>
    std::vector<int16_t> ToInt16Array() const;

    /// \returns tag data as std::vector<uint16_t>
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::vector<uint16_t>
    std::vector<uint16_t> ToUInt16Array() const;

    /// \returns tag data as std::vector<int32_t>
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::vector<int32_t>
    std::vector<int32_t> ToInt32Array() const;

    /// \returns tag data as std::vector<uint32_t>
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::vector<uint32_t>
    std::vector<uint32_t> ToUInt32Array() const;

    /// \returns tag data as std::vector<float>
    /// \throws std::runtime_error if tag does not contain a value of explicit
    ///         type: std::vector<float>
    std::vector<float> ToFloatArray() const;

    /// \}

public:

    /// \name Data Conversion & Validation
    ///

    /// \returns true if tag is null (e.g. default-constructed)
    bool IsNull() const;

    /// \returns true if tag contains a value of type: int8_t
    bool IsInt8() const;

    /// \returns true if tag contains a value of type: uint8_t
    bool IsUInt8() const;

    /// \returns true if tag contains a value of type: int16_t
    bool IsInt16() const;

    /// \returns true if tag contains a value of type: uint16_t
    bool IsUInt16() const;

    /// \returns true if tag contains a value of type: int32_t
    bool IsInt32() const;

    /// \returns true if tag contains a value of type: uint32_t
    bool IsUInt32() const;

    /// \returns true if tag contains a value of type: float
    bool IsFloat() const;

    /// \returns true if tag contains a value of type: std::string
    bool IsString() const;

    /// \returns true if tag contains a value of type: std::string \b AND has a
    ///          TagModifier of TagModifier::HEX_STRING
    bool IsHexString() const;

    /// \returns true if tag contains a value of type: std::vector<int8_t>
    bool IsInt8Array() const;

    /// \returns true if tag contains a value of type: std::vector<uint8_t>
    bool IsUInt8Array() const;

    /// \returns true if tag contains a value of type: std::vector<int16_t>
    bool IsInt16Array() const;

    /// \returns true if tag contains a value of type: std::vector<uint16_t>
    bool IsUInt16Array() const;

    /// \returns true if tag contains a value of type: std::vector<int32_t>
    bool IsInt32Array() const;

    /// \returns true if tag contains a value of type: std::vector<uint32_t>
    bool IsUInt32Array() const;

    /// \returns true if tag contains a value of type: std::vector<float>
    bool IsFloatArray() const;

    /// \returns true if tag contains a value with any signed integer type
    bool IsSignedInt() const;

    /// \returns true if tag contains a value with any unsigned integer type
    bool IsUnsignedInt() const;

    /// \returns true if tag contains a value with any integer type
    bool IsIntegral() const;

    /// \returns true if tag contains a value with any integer or float type
    bool IsNumeric() const;

    /// \returns true if tag contains a vector containing signed integers
    bool IsSignedArray() const;

    /// \returns true if tag contains a vector containing unsigned integers
    bool IsUnsignedArray() const;

    /// \returns true if tag contains a vector containing integers
    bool IsIntegralArray() const;

    /// \returns true if tag contains a vector (integers or floats)
    bool IsArray() const;

    /// \}

public:
    /// \name Type & Modifier Attributes
    /// \{

    /// \returns enum value for current tag data
    TagDataType Type() const;

    /// \returns printable type name for current tag data
    std::string Typename() const;

    /// \returns true if tag data modifier \p m is set
    bool HasModifier(const TagModifier m) const;

    /// \returns current tag data modifier
    TagModifier Modifier() const;

    /// \brief Sets tag data modifier.
    ///
    /// \param[in] m    new modifier value
    ///
    /// \returns reference to this tag
    Tag& Modifier(const TagModifier m);

    /// \}

private :
    // NOTE - keep this synced with TagDataType enum ordering
    using var_t = boost::variant<boost::blank, // <-- default constructor creates variant of this type
                                 int8_t,
                                 uint8_t,
                                 int16_t,
                                 uint16_t,
                                 int32_t,
                                 uint32_t,
                                 float,
                                 std::string,
                                 std::vector<int8_t>,
                                 std::vector<uint8_t>,
                                 std::vector<int16_t>,
                                 std::vector<uint16_t>,
                                 std::vector<int32_t>,
                                 std::vector<uint32_t>,
                                 std::vector<float> >;

    var_t data_;
    TagModifier modifier_ = TagModifier::NONE;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/Tag.inl"

#endif // TAG_H
