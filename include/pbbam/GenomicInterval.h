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

/// This class represents a genomic interval (reference ID, and 0-based coordinates)
class PBBAM_EXPORT GenomicInterval
{
public:

    /// \name Constructors & Related Methods
    ///  \{

    /** Default constructor; yields an empty genomic interval [0,0) with id of -1  */
    GenomicInterval(void);

    /** Constructor for interval from [start, stop) */
    GenomicInterval(const int id,
                    const Position& start,
                    const Position& stop);

    /** Copy constructor */
    GenomicInterval(const GenomicInterval& other);

    ~GenomicInterval(void);

    /// \}

    /// \name Attributes
    /// \{

    /// \returns interval reference id
    /// \sa BamFile::ReferenceId
    inline int Id(void) const;

    /// \returns underlying Interval object
    inline PacBio::BAM::Interval<Position> Interval(void) const;

    /// \returns interval start coordinate
    inline Position Start(void) const;

    /// \returns interval stop coordinate
    inline Position Stop(void) const;

    /// Sets this interval's reference id.
    ///
    /// \param[in] id
    /// \returns reference to this interval
    /// \sa BamFile::ReferenceId
    inline GenomicInterval& Id(const int id);

    /// Sets this underlying Interval
    ///
    /// \param[in] interval
    /// \returns reference to this interval
    inline GenomicInterval& Interval(const PacBio::BAM::Interval<Position>& interval);

    /// Sets this interval's start coordinate.
    ///
    /// \param[in] start
    /// \returns reference to this interval
    inline GenomicInterval& Start(const Position start);

    /// Sets this interval's stop coordinate.
    ///
    /// \param[in] stop
    /// \returns reference to this interval
    inline GenomicInterval& Stop(const Position stop);

    /// \}

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

    /// \name Comparison Operators
    /// \{

    /// \returns true if same id & underlying interval
    inline bool operator==(const GenomicInterval& other) const;

    /// \returns true if either ids or underlying intervals differ
    inline bool operator!=(const GenomicInterval& other) const;

    /// \}

private:
    int id_;
    PacBio::BAM::Interval<Position> interval_;
};

inline GenomicInterval::~GenomicInterval(void) { }

inline int GenomicInterval::Id(void) const
{ return id_; }

inline GenomicInterval& GenomicInterval::Id(const int id)
{ id_ = id; return *this; }

inline PacBio::BAM::Interval<Position> GenomicInterval::Interval(void) const
{ return interval_; }

inline GenomicInterval& GenomicInterval::Interval(const PacBio::BAM::Interval<Position>& interval)
{ interval_ = interval; return *this; }

inline bool GenomicInterval::IsValid(void) const
{
    return id_ >= 0 &&
           interval_.Start() >= 0 &&
           interval_.Stop()  >= 0 &&
           interval_.IsValid();
}

inline Position GenomicInterval::Start(void) const
{ return interval_.Start(); }

inline GenomicInterval& GenomicInterval::Start(const Position start)
{ interval_.Start(start); return *this; }

inline Position GenomicInterval::Stop(void) const
{ return interval_.Stop(); }

inline GenomicInterval& GenomicInterval::Stop(const Position stop)
{ interval_.Stop(stop); return *this; }

inline bool GenomicInterval::operator==(const GenomicInterval& other) const
{ return id_ == other.id_ && interval_ == other.interval_; }

inline bool GenomicInterval::operator!=(const GenomicInterval& other) const
{ return !(*this == other); }

} // namespace BAM
} // namspace PacBio

#endif // GENOMICINTERVAL_H
