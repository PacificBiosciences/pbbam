// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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
#include <string>

namespace PacBio {
namespace BAM {

/// This class represents a genomic interval (reference ID, and 0-based coordinates)
class PBBAM_EXPORT GenomicInterval
{
public:
    typedef typename Interval<uint32_t>::difference_type difference_type;

public:

    /// \name Constructors
    ///  \{

    /** Default constructor; yields an empty genomic interval [0,0) */
    GenomicInterval(void);

    /** Constructor for interval from [start, stop) */
    GenomicInterval(const int id,
                    const int32_t& start,
                    const int32_t& stop);

    /// \}

    /// \name Attributes
    /// \{

    inline int Id(void) const;
    inline GenomicInterval& Id(const int id);

    inline uint32_t Start(void) const;
    inline GenomicInterval& Start(const uint32_t start);

    inline uint32_t Stop(void) const;
    inline GenomicInterval& Stop(const uint32_t stop);

    /// \}

    /// \name Interval Operations
    /// \{

    bool CoveredBy(const GenomicInterval& other) const;
    bool Covers(const GenomicInterval& other) const;
    bool Intersects(const GenomicInterval& other) const;
    difference_type Length(void) const;

    /// \}

    /// \name Comparison Operators
    /// \{

    inline bool operator==(const GenomicInterval& other) const;
    inline bool operator!=(const GenomicInterval& other) const;

    /// \}

private:
    int id_;
    Interval<uint32_t> interval_;
};

inline int GenomicInterval::Id(void) const
{
    return id_;
}

inline GenomicInterval& GenomicInterval::Id(const int id)
{
    id_ = id;
    return *this;
}

inline uint32_t GenomicInterval::Start(void) const
{
    return interval_.Start();
}

inline GenomicInterval& GenomicInterval::Start(const uint32_t start)
{
    interval_.Start(start);
    return *this;
}

inline uint32_t GenomicInterval::Stop(void) const
{
    return interval_.Stop();
}

inline GenomicInterval& GenomicInterval::Stop(const uint32_t stop)
{
    interval_.Stop(stop);
    return *this;
}

inline bool GenomicInterval::operator==(const GenomicInterval& other) const
{
    return id_ == other.id_ && interval_ == other.interval_;
}

inline bool GenomicInterval::operator!=(const GenomicInterval& other) const
{
   return !(*this == other);
}

} // namespace BAM
} // namspace PacBio

#endif // GENOMICINTERVAL_H
