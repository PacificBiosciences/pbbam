// File Description
/// \file Interval.inl
/// \brief Inline implementations for the Interval class.
//
// Author: Derek Barnett

#include "pbbam/Interval.h"

namespace PacBio {
namespace BAM {

template<typename T>
inline Interval<T>::Interval()
    : data_{boost::icl::discrete_interval<T>::right_open(0,0)}
{ }

template<typename T>
inline Interval<T>::Interval(const T val)
    : data_{boost::icl::discrete_interval<T>::right_open(val,val+1)}
{ }

template<typename T>
inline Interval<T>::Interval(const T start, const T stop)
    : data_{boost::icl::discrete_interval<T>::right_open(start,stop)}
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
inline bool Interval<T>::IsValid() const
{ return !boost::icl::is_empty(data_); }

template<typename T>
inline size_t Interval<T>::Length() const
{ return boost::icl::length(data_); }

template<typename T>
inline T Interval<T>::Start() const
{ return data_.lower(); }

template<typename T>
inline Interval<T>& Interval<T>::Start(const T& start)
{
    data_ = boost::icl::discrete_interval<T>::right_open(start, data_.upper());
    return *this;
}

template<typename T>
inline T Interval<T>::Stop() const
{ return data_.upper(); }

template<typename T>
inline Interval<T>& Interval<T>::Stop(const T& stop)
{
    data_ = boost::icl::discrete_interval<T>::right_open(data_.lower(), stop);
    return *this;
}

} // namespace BAM
} // namespace PacBio
