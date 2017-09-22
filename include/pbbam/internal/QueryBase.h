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

#ifndef QUERYBASE_H
#define QUERYBASE_H

#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/DataSet.h"
#include <memory>
#include <vector>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

template<typename T>
class QueryBase;

template<typename T>
class QueryIteratorBase
{
public:
    ~QueryIteratorBase() = default;

    bool operator==(const QueryIteratorBase<T>& other) const;
    bool operator!=(const QueryIteratorBase<T>& other) const;

protected:
    QueryIteratorBase() = default;
    QueryIteratorBase(QueryBase<T>& query);

    void ReadNext();

protected:
    QueryBase<T>* query_ = nullptr;
    T record_;
};

template<typename T>
class QueryIterator : public QueryIteratorBase<T>
{
public:
    QueryIterator() = default;
    QueryIterator(QueryBase<T>& query);

    T& operator*();
    T* operator->();

    QueryIterator<T>& operator++();
    QueryIterator<T> operator++(int);
};

template<typename T>
class QueryConstIterator : public QueryIteratorBase<T>
{
public:
    QueryConstIterator() = default;
    QueryConstIterator(const QueryBase<T>& query);

    const T& operator*() const;
    const T* operator->() const;

    QueryConstIterator<T>& operator++();
    QueryConstIterator<T> operator++(int);
};

template<typename T>
class QueryBase {

public:
    using iterator = QueryIterator<T>;
    using const_iterator = QueryConstIterator<T>;

public:
    virtual ~QueryBase() = default;

public:
    QueryConstIterator<T> begin() const;
    QueryConstIterator<T> cbegin() const;
    QueryIterator<T> begin();

    QueryConstIterator<T> end() const;
    QueryConstIterator<T> cend() const;
    QueryIterator<T> end();

public:
    virtual bool GetNext(T& r) =0;

protected:
    QueryBase() = default;
};

using IQuery = QueryBase<BamRecord>;
using IGroupQuery = QueryBase<std::vector<BamRecord>>;

} // namespace internal
} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/QueryBase.inl"

#endif // QUERYBASE_H
