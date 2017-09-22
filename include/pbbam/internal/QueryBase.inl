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

#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {
namespace internal {

// -------------------
// QueryIteratorBase
// -------------------

template<typename T>
inline QueryIteratorBase<T>::QueryIteratorBase(QueryBase<T>& query)
    : query_(&query)
{ ReadNext(); }

template<typename T> inline
bool QueryIteratorBase<T>::operator==(const QueryIteratorBase<T>& other) const
{ return query_ == other.query_; }

template<typename T> inline
bool QueryIteratorBase<T>::operator!=(const QueryIteratorBase<T>& other) const
{ return !(*this == other); }

// -------------------
// QueryIterator
// -------------------

template<typename T> inline
QueryIterator<T>::QueryIterator(QueryBase<T>& query)
    : QueryIteratorBase<T>(query)
{ }

template<typename T> inline
T& QueryIterator<T>::operator*()
{ return QueryIteratorBase<T>::record_; }

template<typename T> inline
T* QueryIterator<T>::operator->()
{ return &(operator*()); }

template<typename T> inline
QueryIterator<T>& QueryIterator<T>::operator++()
{ QueryIteratorBase<T>::ReadNext(); return *this; }

template<typename T> inline
QueryIterator<T> QueryIterator<T>::operator++(int)
{
    QueryIterator<T> result(*this);
    ++(*this);
    return result;
}

// --------------------
// QueryConstIterator
// --------------------

template<typename T> inline
QueryConstIterator<T>::QueryConstIterator(const QueryBase<T>& query)
    : QueryIteratorBase<T>(const_cast<QueryBase<T>&>(query))
{ }

template<typename T> inline
const T& QueryConstIterator<T>::operator*() const
{ return QueryIteratorBase<T>::record_; }

template<typename T> inline
const T* QueryConstIterator<T>::operator->() const
{ return &(operator*()); }

template<typename T> inline
QueryConstIterator<T>& QueryConstIterator<T>::operator++()
{ QueryIteratorBase<T>::ReadNext(); return *this; }

template<typename T> inline
QueryConstIterator<T> QueryConstIterator<T>::operator++(int)
{
    QueryConstIterator<T> result(*this);
    ++(*this);
    return result;
}

// -----------
// QueryBase
// -----------

template<typename T> inline
QueryConstIterator<T> QueryBase<T>::begin() const
{ return QueryConstIterator<T>(*this); }

template<typename T> inline
QueryIterator<T> QueryBase<T>::begin()
{ return QueryIterator<T>(*this); }

template<typename T> inline
QueryConstIterator<T> QueryBase<T>::cbegin() const
{ return QueryConstIterator<T>(*this); }

template<typename T> inline
QueryConstIterator<T> QueryBase<T>::cend() const
{ return QueryConstIterator<T>(); }

template<typename T> inline
QueryConstIterator<T> QueryBase<T>::end() const
{ return QueryConstIterator<T>(); }

template<typename T> inline
QueryIterator<T> QueryBase<T>::end()
{ return QueryIterator<T>(); }

template<typename T>
inline void QueryIteratorBase<T>::ReadNext()
{
    assert(query_);
    if (!query_->GetNext(record_))
        query_ = nullptr;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio
