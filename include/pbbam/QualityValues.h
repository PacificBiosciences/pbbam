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

// Author: Derek Barnett

#ifndef QUALITYVALUES_H
#define QUALITYVALUES_H

#include "pbbam/QualityValue.h"
#include <algorithm>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The QualityValues class represents a sequence of FASTQ-compatible
/// quality values. See QualityValue documentation for details.
///
class PBBAM_EXPORT QualityValues : public std::vector<QualityValue>
{
public:
    /// Creates a QualityValues collection from a FASTQ-encoded string.
    static QualityValues FromFastq(const std::string& fastq);

public:
    /// \name Constructors & Related Methods
    ///  \{

    QualityValues(void);
    QualityValues(const std::vector<QualityValue>& quals);
    QualityValues(const std::vector<uint8_t>& quals);

    QualityValues(const std::vector<uint8_t>::const_iterator first,
                  const std::vector<uint8_t>::const_iterator last);
    QualityValues(const QualityValues::const_iterator first,
                  const QualityValues::const_iterator last);

    QualityValues(const QualityValues& other);

    QualityValues(std::vector<QualityValue>&& quals);
    QualityValues(QualityValues&& other);

    QualityValues& operator=(const QualityValues& other);
    QualityValues& operator=(const std::vector<QualityValue>& quals);

    QualityValues& operator=(QualityValues&& other);
    QualityValues& operator=(std::vector<QualityValue>&& quals);

    ~QualityValues(void);

    /// \}

public:
    /// \name Iterators
    /// \{

    /// \returns A const_iterator to the beginning of the sequence.
    std::vector<QualityValue>::const_iterator cbegin(void) const;

    /// \returns A const_iterator to the element past the end of the sequence.
    std::vector<QualityValue>::const_iterator cend(void) const;

    /// \returns A const_iterator to the beginning of the sequence.
    std::vector<QualityValue>::const_iterator begin(void) const;

    /// \returns A const_iterator to the element past the end of the sequence.
    std::vector<QualityValue>::const_iterator end(void) const;

    /// \returns An iterator to the beginning of the sequence.
    std::vector<QualityValue>::iterator begin(void);

    /// \returns An iterator to the element past the end of the sequence.
    std::vector<QualityValue>::iterator end(void);

    /// \}

public:
    /// \returns the FASTQ-encoded string for this collection
    std::string Fastq(void) const;
};

inline QualityValues::QualityValues(void)
    : std::vector<QualityValue>()
{ }

inline QualityValues::QualityValues(const std::vector<QualityValue>& quals)
    : std::vector<QualityValue>(quals)
{ }

inline QualityValues::QualityValues(const std::vector<uint8_t>& quals)
    : std::vector<QualityValue>()
{
    resize(quals.size());
    std::copy(quals.cbegin(), quals.cend(), begin());
}

inline QualityValues::QualityValues(const std::vector<uint8_t>::const_iterator first,
                                    const std::vector<uint8_t>::const_iterator last)
    : std::vector<QualityValue>(first, last)
{ }

inline QualityValues::QualityValues(const QualityValues::const_iterator first,
                                    const QualityValues::const_iterator last)
    : std::vector<QualityValue>()
{
    assign(first, last);
}

inline QualityValues::QualityValues(const QualityValues& other)
    : std::vector<QualityValue>(other)
{ }

inline QualityValues::QualityValues(std::vector<QualityValue>&& quals)
    : std::vector<QualityValue>(std::move(quals))
{ }

inline QualityValues::QualityValues(QualityValues&& other)
    : std::vector<QualityValue>(std::move(other))
{ }

inline QualityValues& QualityValues::operator=(const QualityValues& other)
{ std::vector<QualityValue>::operator=(other); return *this; }

inline QualityValues& QualityValues::operator=(const std::vector<QualityValue>& quals)
{ std::vector<QualityValue>::operator=(quals); return *this; }

inline QualityValues& QualityValues::operator=(QualityValues&& other)
{ std::vector<QualityValue>::operator=(std::move(other)); return *this; }

inline QualityValues& QualityValues::operator=(std::vector<QualityValue>&& quals)
{ std::vector<QualityValue>::operator=(std::move(quals)); return *this; }

inline QualityValues::~QualityValues(void) { }

inline std::vector<QualityValue>::const_iterator QualityValues::cbegin(void) const
{ return std::vector<QualityValue>::cbegin(); }

inline std::vector<QualityValue>::const_iterator QualityValues::cend(void) const
{ return std::vector<QualityValue>::cend(); }

inline std::vector<QualityValue>::const_iterator QualityValues::begin(void) const
{ return std::vector<QualityValue>::begin(); }

inline std::vector<QualityValue>::const_iterator QualityValues::end(void) const
{ return std::vector<QualityValue>::end(); }

inline std::vector<QualityValue>::iterator QualityValues::begin(void)
{ return std::vector<QualityValue>::begin(); }

inline std::vector<QualityValue>::iterator QualityValues::end(void)
{ return std::vector<QualityValue>::end(); }

inline QualityValues QualityValues::FromFastq(const std::string& fastq)
{
    QualityValues result;
    result.resize(fastq.size());
    std::transform(fastq.cbegin(), fastq.cend(), result.begin(), QualityValue::FromFastq);
    return result;
}

inline std::string QualityValues::Fastq(void) const
{
    std::string result;
    result.reserve(size());
    auto iter = cbegin();
    const auto end = cend();
    for (; iter != end; ++iter)
        result.push_back((*iter).Fastq());
    return result;
}

inline bool operator==(const QualityValues& lhs, const std::string& rhs)
{ return lhs == QualityValues::FromFastq(rhs); }

inline bool operator==(const std::string& lhs, const QualityValues& rhs)
{ return rhs == lhs; }

inline bool operator!=(const QualityValues& lhs, const std::string& rhs)
{ return !(lhs == rhs); }

inline bool operator!=(const std::string& lhs, const QualityValues& rhs)
{ return !(lhs == rhs); }

} // namespace BAM
} // namespace PacBio

#endif // QUALITYVALUES_H
