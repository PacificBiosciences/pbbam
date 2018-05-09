// File Description
/// \file Interval.h
/// \brief Defines the Interval class.
//
// Author: Derek Barnett

#ifndef INTERVAL_H
#define INTERVAL_H

#include <cstddef>
#include <string>
#include "pbbam/Config.h"

#define BOOST_ICL_USE_STATIC_BOUNDED_INTERVALS
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/interval_traits.hpp>

namespace PacBio {
namespace BAM {

/// \brief Represents a half-open (right-open) interval [start, stop)
///
/// \note This class is agnostic whether the values are 0-based or 1-based.
///       Client code should primarily work with GenomicInterval, which does
///       enforce this distinction.
///
template <typename T>
class Interval
{
public:
    using interval_type = boost::icl::discrete_interval<T>;

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty interval [0,0)
    Interval();

    /// \brief Creates a 'singleton' interval [val,val+1)
    Interval(const T val);

    /// brief Creates an interval from [start, stop) */
    Interval(const T start, const T stop);

    Interval(const Interval<T>&) = default;
    Interval(Interval&&) = default;
    Interval& operator=(const Interval<T>&) = default;
    Interval& operator=(Interval<T>&&) = default;
    ~Interval() = default;

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    /// \returns true if both intervals share the same endpoints
    bool operator==(const Interval<T>& other) const;

    /// \returns true if either interval's endpoints differ
    bool operator!=(const Interval<T>& other) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \returns interval's start coordinate
    T Start() const;

    /// Sets this interval's start coordinate.
    ///
    /// \param[in] start
    /// \returns reference to this interval
    ///
    Interval<T>& Start(const T& start);

    /// \returns interval's stop coordinate
    T Stop() const;

    /// Sets this interval's stop coordinate.
    ///
    /// \param[in] stop
    /// \returns reference to this interval
    ///
    Interval<T>& Stop(const T& stop);

    /// \}

public:
    /// \name Interval Operations

    /// \returns true if this interval is fully covered by (or contained in) \p other
    bool CoveredBy(const Interval<T>& other) const;

    //// \returns true if this interval covers (or contains) \p other
    bool Covers(const Interval<T>& other) const;

    /// \returns true if intervals interset
    bool Intersects(const Interval<T>& other) const;

    /// \returns true if interval is valid (e.g. start < stop)
    bool IsValid() const;

    /// \returns interval length
    size_t Length() const;

    /// \}

private:
    interval_type data_;
};

}  // namespace BAM
}  // namspace PacBio

#include "pbbam/internal/Interval.inl"

#endif  // GENOMICINTERVAL_H
