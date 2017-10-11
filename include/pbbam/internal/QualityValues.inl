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
/// \file QualityValues.inl
/// \brief Inline implementations for the QualityValues class.
//
// Author: Derek Barnett

#include "pbbam/QualityValues.h"
#include <algorithm>

namespace PacBio {
namespace BAM {

inline QualityValues::QualityValues()
    : std::vector<QualityValue>()
{ }

inline QualityValues::QualityValues(const std::string& fastqString)
    : std::vector<QualityValue>()
{
    resize(fastqString.size());
    std::transform(fastqString.cbegin(), fastqString.cend(),
                   begin(), QualityValue::FromFastq);
}

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

inline QualityValues::QualityValues(std::vector<QualityValue>&& quals)
    : std::vector<QualityValue>(std::move(quals))
{ }

inline QualityValues& QualityValues::operator=(const std::vector<QualityValue>& quals)
{ std::vector<QualityValue>::operator=(quals); return *this; }

inline QualityValues& QualityValues::operator=(std::vector<QualityValue>&& quals)
{ std::vector<QualityValue>::operator=(std::move(quals)); return *this; }

inline std::vector<QualityValue>::const_iterator QualityValues::cbegin() const
{ return std::vector<QualityValue>::cbegin(); }

inline std::vector<QualityValue>::const_iterator QualityValues::cend() const
{ return std::vector<QualityValue>::cend(); }

inline std::vector<QualityValue>::const_iterator QualityValues::begin() const
{ return std::vector<QualityValue>::begin(); }

inline std::vector<QualityValue>::const_iterator QualityValues::end() const
{ return std::vector<QualityValue>::end(); }

inline std::vector<QualityValue>::iterator QualityValues::begin()
{ return std::vector<QualityValue>::begin(); }

inline std::vector<QualityValue>::iterator QualityValues::end()
{ return std::vector<QualityValue>::end(); }

inline QualityValues QualityValues::FromFastq(const std::string& fastq)
{ return QualityValues(fastq); }

inline std::string QualityValues::Fastq() const
{
    std::string result;
    result.reserve(size());
    auto iter = cbegin();
    const auto end = cend();
    for (; iter != end; ++iter)
        result.push_back((*iter).Fastq());
    return result;
}

inline bool QualityValues::operator==(const std::string& fastq) const
{ return *this == QualityValues(fastq); }

inline bool QualityValues::operator!=(const std::string& fastq) const
{ return *this != QualityValues(fastq); }

} // namespace BAM
} // namespace PacBio
