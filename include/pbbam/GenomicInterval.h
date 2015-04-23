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

#ifndef GENOMICINTERVAL_H
#define GENOMICINTERVAL_H

#include "pbbam/Config.h"
#include "pbbam/Interval.h"
#include "pbbam/Position.h"
#include <string>

namespace PacBio {
namespace BAM {

/// This class represents a genomic interval (reference name, and 0-based coordinates)
class PBBAM_EXPORT GenomicInterval
{
public:
    /// \name Constructors & Related Methods
    ///  \{

    /** Default constructor; yields an empty genomic interval */
    GenomicInterval(void);

    /** Constructor for interval on named sequence with range: [start, stop) */
    GenomicInterval(const std::string& name,
                    const Position& start,
                    const Position& stop);

    /** Copy constructor */
    GenomicInterval(const GenomicInterval& other);

    ~GenomicInterval(void);

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \returns interval reference name
    std::string Name(void) const;

    /// \returns underlying Interval object
    PacBio::BAM::Interval<Position> Interval(void) const;

    /// \returns interval start coordinate
    Position Start(void) const;

    /// \returns interval stop coordinate
    Position Stop(void) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// Sets this interval's reference name.
    ///
    /// \param[in] name
    /// \returns reference to this interval
    GenomicInterval& Name(const std::string& name);

    /// Sets this underlying Interval
    ///
    /// \param[in] interval
    /// \returns reference to this interval
    GenomicInterval& Interval(const PacBio::BAM::Interval<Position>& interval);

    /// Sets this interval's start coordinate.
    ///
    /// \param[in] start
    /// \returns reference to this interval
    GenomicInterval& Start(const Position start);

    /// Sets this interval's stop coordinate.
    ///
    /// \param[in] stop
    /// \returns reference to this interval
    GenomicInterval& Stop(const Position stop);

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

    /// \returns true if underlying Interval::IsValid(), and id/endpoints are non-negative.
    bool IsValid(void) const;

    /// \returns length of underlying
    size_t Length(void) const;

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    /// \returns true if same id & underlying interval
    bool operator==(const GenomicInterval& other) const;

    /// \returns true if either ids or underlying intervals differ
    bool operator!=(const GenomicInterval& other) const;

    /// \}

private:
    std::string name_;
    PacBio::BAM::Interval<Position> interval_;
};

inline GenomicInterval::~GenomicInterval(void) { }

inline std::string GenomicInterval::Name(void) const
{ return name_; }

inline GenomicInterval& GenomicInterval::Name(const std::string& name)
{ name_ = name; return *this; }

inline PacBio::BAM::Interval<Position> GenomicInterval::Interval(void) const
{ return interval_; }

inline GenomicInterval& GenomicInterval::Interval(const PacBio::BAM::Interval<Position>& interval)
{ interval_ = interval; return *this; }

inline bool GenomicInterval::IsValid(void) const
{
    return !name_.empty() &&
           interval_.Start() >= 0 &&
           interval_.Stop()  >= 0 &&
           interval_.IsValid();
}

inline size_t GenomicInterval::Length(void) const
{ return interval_.Length(); }

inline Position GenomicInterval::Start(void) const
{ return interval_.Start(); }

inline GenomicInterval& GenomicInterval::Start(const Position start)
{ interval_.Start(start); return *this; }

inline Position GenomicInterval::Stop(void) const
{ return interval_.Stop(); }

inline GenomicInterval& GenomicInterval::Stop(const Position stop)
{ interval_.Stop(stop); return *this; }

inline bool GenomicInterval::operator==(const GenomicInterval& other) const
{ return name_ == other.name_ && interval_ == other.interval_; }

inline bool GenomicInterval::operator!=(const GenomicInterval& other) const
{ return !(*this == other); }

} // namespace BAM
} // namspace PacBio

#endif // GENOMICINTERVAL_H
