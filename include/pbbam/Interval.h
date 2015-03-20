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

#ifndef INTERVAL_H
#define INTERVAL_H

#include "pbbam/Config.h"
#include <string>

#define BOOST_ICL_USE_STATIC_BOUNDED_INTERVALS
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/interval_traits.hpp>

namespace PacBio {
namespace BAM {

/// \brief Utility class for working with half-open (right-open) intervals. [start, stop)
///
/// \note This class is agnostic whether the values are 0-based or 1-based.
/// \todo Should it be? Should that go here or "higher up"?
///
template<typename T>
class Interval
{
public:
    typedef boost::icl::discrete_interval<T> interval_type;
    typedef typename boost::icl::interval_traits<interval_type> interval_traits;
    typedef typename boost::icl::difference_type_of<interval_traits>::type difference_type;

public:

    /// \name Constructors
    /// \{

    /** Default constructor; yields an empty interval [0,0) */
    inline Interval(void);

    /** Constructor for a singleton interval [val,val+1) */
    inline Interval(const T& val);

    /** Constructor for interval from [start, stop) */
    inline Interval(const T& start, const T& stop);

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \returns interval start coordinate
    inline T Start(void) const;

    /// Sets this interval's start coordinate.
    ///
    /// \param[in] start
    /// \returns reference to this interval
    inline Interval<T>& Start(const T& start);

    /// \returns interval stop coordinate
    inline T Stop(void) const;

    /// Sets this interval's stop coordinate.
    ///
    /// \param[in] stop
    /// \returns reference to this interval
    inline Interval<T>& Stop(const T& stop);

    /// \}

    /// \name Interval Operations

    /// \returns true if this interval is fully covered by (or contained in) \p other
    inline bool CoveredBy(const Interval<T>& other) const;

    //// \returns true if this interval covers (or contains) \p other
    inline bool Covers(const Interval<T>& other) const;

    /// \returns true if intervals interset
    inline bool Intersects(const Interval<T>& other) const;

    /// \returns true if interval is valid (e.g. start < stop)
    inline bool IsValid(void) const;

    /// \returns interval length
    inline difference_type Length(void) const;

    /// \}

    /// \name Comparison Operators
    /// \{

    /// \returns true if both intervals share the same endpoints
    inline bool operator==(const Interval<T>& other) const;

    /// \returns true if either interval's endpoints differ
    inline bool operator!=(const Interval<T>& other) const;

    /// \}

private:
    interval_type data_;
};

template<typename T>
Interval<T>::Interval(void)
    : data_(boost::icl::discrete_interval<T>::right_open(0,0))
{ }

template<typename T>
Interval<T>::Interval(const T& val)
    : data_(boost::icl::discrete_interval<T>::right_open(val,val+1))
{ }

template<typename T>
Interval<T>::Interval(const T& start, const T& stop)
    : data_(boost::icl::discrete_interval<T>::right_open(start,stop))
{ }

template<typename T>
inline bool Interval<T>::operator==(const Interval<T>& other) const
{ return data_ == other.data_; }

template<typename T>
inline bool Interval<T>::operator!=(const Interval<T>& other) const
{ return !(data_ == other.data_); }

template<typename T>
inline bool Interval<T>::CoveredBy(const Interval<T>& other) const
{ return boost::icl::within(data_, other.data_); }

template<typename T>
inline bool Interval<T>::Covers(const Interval<T>& other) const
{ return boost::icl::contains(data_, other.data_); }

template<typename T>
inline bool Interval<T>::Intersects(const Interval<T>& other) const
{ return boost::icl::intersects(data_, other.data_); }

template<typename T>
inline bool Interval<T>::IsValid(void) const
{ return !boost::icl::is_empty(data_); }

template<typename T>
inline typename Interval<T>::difference_type Interval<T>::Length(void) const
{ return boost::icl::length(data_); }

template<typename T>
inline T Interval<T>::Start(void) const
{ return data_.lower(); }

template<typename T>
inline Interval<T>& Interval<T>::Start(const T& start)
{
    data_ = boost::icl::discrete_interval<T>::right_open(start, data_.upper());
    return *this;
}

template<typename T>
inline T Interval<T>::Stop(void) const
{ return data_.upper(); }

template<typename T>
inline Interval<T>& Interval<T>::Stop(const T& stop)
{
    data_ = boost::icl::discrete_interval<T>::right_open(data_.lower(), stop);
    return *this;
}

} // namespace BAM
} // namspace PacBio

#endif // GENOMICINTERVAL_H
