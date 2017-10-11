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

#include "pbbam/internal/DataSetListElement.h"
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

// --------------------
// DataSetListElement
// --------------------

template<class T>
inline DataSetListElement<T>::DataSetListElement(const std::string& label,
                                                 const XsdType& xsd)
    : DataSetElement(label, xsd)
{ }

template<class T>
inline const T& DataSetListElement<T>::operator[](size_t index) const
{ return static_cast<const T&>(children_.at(index)); }

template<class T>
inline T& DataSetListElement<T>::operator[](size_t index)
{ return static_cast<T&>(children_.at(index)); }

template<class T>
inline size_t DataSetListElement<T>::Size() const
{ return NumChildren(); }

template<class T>
inline DataSetListIterator<T> DataSetListElement<T>::begin()
{ return DataSetListIterator<T>(this, 0); }

template<class T>
inline DataSetListConstIterator<T> DataSetListElement<T>::begin() const
{ return DataSetListConstIterator<T>(this, 0); }

template<class T>
inline DataSetListConstIterator<T> DataSetListElement<T>::cbegin() const
{ return DataSetListConstIterator<T>(this, 0); }

template<class T>
inline DataSetListIterator<T> DataSetListElement<T>::end()
{ return DataSetListIterator<T>(this, NumChildren()); }

template<class T>
inline DataSetListConstIterator<T> DataSetListElement<T>::end() const
{ return DataSetListConstIterator<T>(this, NumChildren()); }

template<class T>
inline DataSetListConstIterator<T>DataSetListElement<T>::cend() const
{ return DataSetListConstIterator<T>(this, NumChildren()); }

// -------------------------
// DataSetListIteratorBase
// -------------------------

template<class T>
inline bool DataSetListIteratorBase<T>::operator==(const DataSetListIteratorBase<T>& other) const
{ return parent_ == other.parent_ &&
         index_ == other.index_;
}

template<class T>
inline bool DataSetListIteratorBase<T>::operator!=(const DataSetListIteratorBase<T>& other) const
{ return !(*this == other); }

template<class T>
inline DataSetListIteratorBase<T>::DataSetListIteratorBase(const DataSetListElement<T>* parent, size_t i)
    : parent_(parent)
    , index_(i)
{ }

template<class T>
inline void DataSetListIteratorBase<T>::ReadNext()
{
    if (index_ >= parent_->NumChildren()) {
        parent_ = nullptr;
        return;
    }
    ++index_;
}

// ---------------------
// DataSetListIterator
// ---------------------

template<class T>
inline DataSetListIterator<T>::DataSetListIterator(const DataSetListElement<T>* parent, size_t i)
    : DataSetListIteratorBase<T>(parent, i)
{ }

template<class T>
inline T& DataSetListIterator<T>::operator*()
{ return DataSetListIteratorBase<T>::parent_->template Child<T>(DataSetListIteratorBase<T>::index_); }

template<class T>
inline T* DataSetListIterator<T>::operator->()
{ return &(operator*()); }

template<class T>
inline DataSetListIterator<T>& DataSetListIterator<T>::operator++()
{ DataSetListIteratorBase<T>::ReadNext(); return *this; }

template<class T>
inline DataSetListIterator<T> DataSetListIterator<T>::operator++(int)
{
    DataSetListIterator<T> result(*this);
    ++(*this);
    return result;
}

// --------------------------
// DataSetListConstIterator
// --------------------------

template<class T>
inline DataSetListConstIterator<T>::DataSetListConstIterator(const DataSetListElement<T>* parent, size_t i)
    : DataSetListIteratorBase<T>(parent, i)
{ }

template<class T>
inline const T& DataSetListConstIterator<T>::operator*() const
{ return DataSetListIteratorBase<T>::parent_->template Child<T>(DataSetListIteratorBase<T>::index_); }

template<class T>
inline const T* DataSetListConstIterator<T>::operator->() const
{ return &(operator*()); }

template<class T>
inline DataSetListConstIterator<T>& DataSetListConstIterator<T>::operator++()
{ DataSetListIteratorBase<T>::ReadNext(); return *this; }

template<class T>
inline DataSetListConstIterator<T> DataSetListConstIterator<T>::operator++(int)
{
    DataSetListConstIterator<T> result(*this);
    ++(*this);
    return result;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio
