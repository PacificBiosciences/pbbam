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
/// \file QualityValue.h
/// \brief Defines the QualityValue class.
//
// Author: Derek Barnett

#ifndef QUALITYVALUE_H
#define QUALITYVALUE_H

#include "pbbam/Config.h"
#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The QualityValue class represents a FASTQ-compatible quality value.
///
/// Integers are clamped to [0, 93] (corresponding to ASCII printable chars
/// [!-~]).
///
/// Use QualityValue::FromFastq for constructing entries from FASTQ encoding
/// characters. Otherwise, the resulting QualityValue will be interpreted using
/// the character's numeric value (ignoring the FASTQ offset of 33).
///
class PBBAM_EXPORT QualityValue
{
public:
    static const uint8_t MAX;

public:
    /// \name Conversion Methods
    /// \{

    /// \brief Creates a QualityValue from a FASTQ-encoding character.
    ///
    /// \param[in] c    FASTQ character
    /// \returns quality value representing (c - 33)
    ///
    static QualityValue FromFastq(const char c);

    /// \}

public:
    /// \name Constructors & Related Methods
    ///  \{

    /// \brief Creates a QualityValue with specified value.
    ///
    /// \param[in] value    quality value
    ///
    QualityValue(const uint8_t value = 0);

    QualityValue(const QualityValue& other) = default;
    QualityValue(QualityValue&& other) = default;
    QualityValue& operator=(const QualityValue& other) = default;
    QualityValue& operator=(QualityValue&& other) = default;
    ~QualityValue() = default;

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// \returns the FASTQ-encoding char for this QualityValue
    char Fastq() const;

    /// \returns the integer value of this QualityValue
    operator uint8_t() const;

    /// \}

private:
    uint8_t value_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/QualityValue.inl"

#endif // QUALITYVALUE_H
