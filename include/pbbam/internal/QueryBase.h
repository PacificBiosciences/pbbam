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
    virtual ~QueryIteratorBase(void);

    bool operator==(const QueryIteratorBase<T>& other) const;
    bool operator!=(const QueryIteratorBase<T>& other) const;

protected:
    QueryIteratorBase(void);
    QueryIteratorBase(QueryBase<T>& query);

    void ReadNext(void);

protected:
    QueryBase<T>* query_;
    T record_;
};

template<typename T>
class QueryIterator : public QueryIteratorBase<T>
{
public:
    QueryIterator(void);
    QueryIterator(QueryBase<T>& query);

    T& operator*(void);
    T* operator->(void);

    QueryIterator<T>& operator++(void);
    QueryIterator<T> operator++(int);
};

template<typename T>
class QueryConstIterator : public QueryIteratorBase<T>
{
public:
    QueryConstIterator(void);
    QueryConstIterator(const QueryBase<T>& query);

    const T& operator*(void) const;
    const T* operator->(void) const;

    QueryConstIterator<T>& operator++(void);
    QueryConstIterator<T> operator++(int);
};

template<typename T>
class QueryBase {

public:
    typedef QueryIterator<T>      iterator;
    typedef QueryConstIterator<T> const_iterator;

public:
    virtual ~QueryBase(void);

public:
    QueryConstIterator<T> begin(void) const;
    QueryConstIterator<T> cbegin(void) const;
    QueryIterator<T> begin(void);

    QueryConstIterator<T> end(void) const;
    QueryConstIterator<T> cend(void) const;
    QueryIterator<T> end(void);

public:
    virtual bool GetNext(T& r) =0;

protected:
    QueryBase(void);
};

typedef QueryBase<BamRecord>               IQuery;
typedef QueryBase<std::vector<BamRecord> > IGroupQuery;

} // namespace internal
} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/QueryBase.inl"

#endif // QUERYBASE_H
