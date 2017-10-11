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
/// \file GenomicInterval.h
/// \brief Defines the GenomicInterval class.
//
// Author: Derek Barnett

#ifndef GENOMICINTERVAL_H
#define GENOMICINTERVAL_H

#include "pbbam/Config.h"
#include "pbbam/Interval.h"
#include "pbbam/Position.h"
#include <cstddef>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The GenomicInterval class represents a genomic interval (reference
///        name and 0-based coordinates).
///
class PBBAM_EXPORT GenomicInterval
{
public:
    /// \name Constructors & Related Methods
    ///  \{

    /// \brief Creates an empty genomic interval
    GenomicInterval() = default;

    /// \brief Creates a genomic interval on sequence with \p name, using range:
    ///       [\p start, \p stop)
    GenomicInterval(std::string name,
                    Position start,
                    Position stop);

    /// \brief Creates a genomic interval, using REGION string
    ///
    /// "<ref>:<start>-<stop>" ("chr8:200-600")
    ///
    /// \note The htslib/samtools REGION string expects start positions to be
    ///       1-based. However, throughout pbbam (including the rest of this
    ///       class), we stick to 0-based start coordinates. Thus, while the
    ///       syntax matches that of samtools, we are using a 0-based start
    ///       coordinate here.
    ///
    GenomicInterval(const std::string& zeroBasedRegionString);

    GenomicInterval(const GenomicInterval& other) = default;
    GenomicInterval(GenomicInterval&& other) = default;
    GenomicInterval& operator=(const GenomicInterval& other) = default;
    GenomicInterval& operator=(GenomicInterval&& other) = default;
    ~GenomicInterval() = default;

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    /// \returns true if same id & underlying interval
    bool operator==(const GenomicInterval& other) const;

    /// \returns true if either ids or underlying intervals differ
    bool operator!=(const GenomicInterval& other) const;

    /// \}

public:
    /// \name Interval Operations
    /// \{

    /// \returns true if same id and underlying Interval::CoveredBy() other.
    bool CoveredBy(const GenomicInterval& other) const;

    /// \returns true if same id and underlying Interval::Covers() other.
    bool Covers(const GenomicInterval& other) const;

    /// \returns true if same id and underlying Interval::Intersects() other.
    bool Intersects(const GenomicInterval& other) const;

    /// \returns true if underlying Interval::IsValid(), and id/endpoints are
    ///          non-negative.
    ///
    bool IsValid() const;

    /// \returns length of underlying
    size_t Length() const;

    /// \}


public:
    /// \name Attributes
    /// \{

    /// \returns interval reference name
    std::string Name() const;

    /// \returns underlying Interval object
    PacBio::BAM::Interval<Position> Interval() const;

    /// \returns interval start coordinate
    Position Start() const;

    /// \returns interval stop coordinate
    Position Stop() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// Sets this interval's reference name.
    ///
    /// \param[in] name
    /// \returns reference to this interval
    ///
    GenomicInterval& Name(const std::string& name);

    /// Sets this underlying Interval
    ///
    /// \param[in] interval
    /// \returns reference to this interval
    ///
    GenomicInterval& Interval(const PacBio::BAM::Interval<Position>& interval);

    /// Sets this interval's start coordinate.
    ///
    /// \param[in] start
    /// \returns reference to this interval
    ///
    GenomicInterval& Start(const Position start);

    /// Sets this interval's stop coordinate.
    ///
    /// \param[in] stop
    /// \returns reference to this interval
    ///
    GenomicInterval& Stop(const Position stop);

    /// \}

private:
    std::string name_;
    PacBio::BAM::Interval<Position> interval_;
};

} // namespace BAM
} // namspace PacBio

#include "pbbam/internal/GenomicInterval.inl"

#endif // GENOMICINTERVAL_H
