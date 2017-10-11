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
/// \file QualityValues.h
/// \brief Defines the QualityValues class.
//
// Author: Derek Barnett

#ifndef QUALITYVALUES_H
#define QUALITYVALUES_H

#include "pbbam/QualityValue.h"
#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The QualityValues class represents a sequence of FASTQ-compatible
///        quality values. See QualityValue documentation for more details.
///
class PBBAM_EXPORT QualityValues : public std::vector<QualityValue>
{
public:
    /// \brief Creates a QualityValues object from a FASTQ-encoded string.
    ///
    /// \param[in] fastq    FASTQ-encoded string
    /// \returns corresponding QualityValues object
    ///
    static QualityValues FromFastq(const std::string& fastq);

public:
    /// \name Constructors & Related Methods
    ///  \{

    /// \brief Default constructor - creates an empty QualityValues object.
    QualityValues();

    /// \brief Creates a QualityValues object from a FASTQ-encoded string.
    ///
    /// \param[in] fastqString  FASTQ-encoded string
    ///
    explicit QualityValues(const std::string& fastqString);

    /// \brief Creates a QualityValues object from a vector of QualityValue
    ///        elements.
    ///
    /// \param[in] quals    vector of QualityValue elements
    ///
    explicit QualityValues(const std::vector<QualityValue>& quals);

    /// \brief Creates a QualityValues object from a vector of QualityValue
    ///        elements.
    ///
    /// \param[in] quals    vector of QualityValue elements
    ///
    QualityValues(std::vector<QualityValue>&& quals);

    /// \brief Creates a QualityValues object from a vector of (numeric) quality
    ///        values.
    ///
    /// \param[in] quals    vector of quality value numbers
    ///
    explicit QualityValues(const std::vector<uint8_t>& quals);

    /// \brief Creates a QualityValues object from the contents of the range:
    ///        [first, last)
    ///
    /// \param[in] first    input iterator, whose element is a numeric quality
    /// \param[in] last     input iterator, whose element is a numeric quality
    ///
    QualityValues(const std::vector<uint8_t>::const_iterator first,
                  const std::vector<uint8_t>::const_iterator last);

    /// \brief Creates a QualityValues object from the contents of the range:
    ///        [first, last)
    ///
    /// \param[in] first    input iterator, whose element is a QualityValue
    /// \param[in] last     input iterator, whose element is a QualityValue
    ///
    QualityValues(const QualityValues::const_iterator first,
                  const QualityValues::const_iterator last);

    /// \brief Copy constructor
    QualityValues(const QualityValues& other) = default;

    /// \brief Move constructor
    QualityValues(QualityValues&& other) = default;

    /// \brief Copy assignment operator
    ///
    /// \param[in] other    QualityValues object
    ///
    QualityValues& operator=(const QualityValues& other) = default;

    /// \brief Move assignment operator
    ///
    /// \param[in] other    QualityValues object
    ///
    QualityValues& operator=(QualityValues&& other) = default;

    /// \brief Copy assignment operator
    ///
    /// \param[in] quals    vector of QualityValue elements
    ///
    QualityValues& operator=(const std::vector<QualityValue>& quals);

    /// \brief Move assignment operator
    ///
    /// \param[in] quals    vector of QualityValue elements
    ///
    QualityValues& operator=(std::vector<QualityValue>&& quals);

    /// \brief Destructor
    ~QualityValues() = default;

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    bool operator==(const std::string& other) const;
    bool operator!=(const std::string& other) const;

    /// \}

public:
    /// \name Iterators
    /// \{

    /// \returns a const_iterator to the beginning of the sequence
    std::vector<QualityValue>::const_iterator cbegin() const;

    /// \returns a const_iterator to the element following the last element
    std::vector<QualityValue>::const_iterator cend() const;

    /// \returns a const_iterator to the beginning of the sequence
    std::vector<QualityValue>::const_iterator begin() const;

    /// \returns a const_iterator to the element following the last element
    std::vector<QualityValue>::const_iterator end() const;

    /// \returns an iterator to the beginning of the sequence
    std::vector<QualityValue>::iterator begin();

    /// \returns an iterator to the element following the last element
    std::vector<QualityValue>::iterator end();

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// \returns the FASTQ-encoded string for this sequence of quality values
    std::string Fastq() const;

    /// \}
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/QualityValues.inl"

#endif // QUALITYVALUES_H
