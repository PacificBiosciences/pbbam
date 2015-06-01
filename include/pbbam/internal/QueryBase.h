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

#ifndef QUERYBASE2_H
#define QUERYBASE2_H

#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/dataset/DataSet.h"
#include "pbbam/internal/FilterEngine.h"
#include "pbbam/internal/IBamFileIterator.h"
#include "pbbam/internal/IMergeStrategy.h"
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
protected:
    QueryIteratorBase(void);
    QueryIteratorBase(QueryBase<T>& query);

public:
    virtual ~QueryIteratorBase(void) { }

protected:
    void ReadNext(void);

public:
    bool operator==(const QueryIteratorBase<T>& other) const
    { return query_ == other.query_; }

    bool operator!=(const QueryIteratorBase<T>& other) const
    { return !(*this == other); }

protected:
    QueryBase<T>* query_;
    T record_;
};

template<typename T>
class QueryIterator : public QueryIteratorBase<T>
{
public:
    QueryIterator(void) : QueryIteratorBase<T>() { }
    QueryIterator(QueryBase<T>& query)
        : QueryIteratorBase<T>(query)
    { }

    T& operator*(void) { return QueryIteratorBase<T>::record_; }
    T* operator->(void) { return &(operator*()); }

    QueryIterator<T>& operator++(void)
    { QueryIteratorBase<T>::ReadNext(); return *this; }

    QueryIterator<T> operator++(int)
    {
        QueryIterator<T> result(*this);
        ++(*this);
        return result;
    }
};

template<typename T>
class QueryConstIterator : public QueryIteratorBase<T>
{
public:
    QueryConstIterator(void) : QueryIteratorBase<T>() { }
    QueryConstIterator(const QueryBase<T>& query)
        : QueryIteratorBase<T>(const_cast<QueryBase<T>&>(query))
    { }

    const T& operator*(void) const { return QueryIteratorBase<T>::record_; }
    const T* operator->(void) const { return &(operator*()); }

    QueryConstIterator<T>& operator++(void)
    { QueryIteratorBase<T>::ReadNext(); return *this; }

    QueryConstIterator<T> operator++(int)
    {
        QueryConstIterator<T> result(*this);
        ++(*this);
        return result;
    }
};

template<typename T>
class QueryBase {

public:
    using FileIterPtr = typename IBamFileIteratorBase<T>::Ptr;

protected:
    QueryBase(const DataSet& dataset);
public:
    virtual ~QueryBase(void) { }

    QueryConstIterator<T> begin(void) const  { return QueryConstIterator<T>(*this); }
    QueryConstIterator<T> cbegin(void) const { return QueryConstIterator<T>(*this); }
    QueryIterator<T> begin(void) { return QueryIterator<T>(*this); }

    QueryConstIterator<T> end(void) const { return QueryConstIterator<T>(); }
    QueryConstIterator<T> cend(void) const { return QueryConstIterator<T>(); }
    QueryIterator<T> end(void) { return QueryIterator<T>(); }

public:
    bool GetNext(T& r);

    std::vector<BamFile> GetBamFiles(void) const
    { return dataset_.BamFiles(); }

public:
    std::vector<FileIterPtr> CreateIterators(void)
    {
        const std::vector<BamFile>& bamFiles = dataset_.BamFiles();
        std::vector<FileIterPtr> result;
        result.reserve(bamFiles.size());
        for ( const BamFile& bamFile : bamFiles)
            result.push_back(CreateIterator(bamFile));
        return result;
    }

protected:
    virtual FileIterPtr CreateIterator(const BamFile& bamFile) = 0;

protected:
    DataSet dataset_;
    std::unique_ptr<IMergeStrategyBase<T> > mergeStrategy_;
    FilterEngine filterEngine_;
};

typedef QueryBase<BamRecord>               IQuery;
typedef QueryBase<std::vector<BamRecord> > IGroupQuery;

template<typename T>
inline QueryIteratorBase<T>::QueryIteratorBase(void)
    : query_(nullptr)
{ }

template<typename T>
inline QueryIteratorBase<T>::QueryIteratorBase(QueryBase<T> &query)
    : query_(&query)
{ ReadNext(); }

template<typename T>
inline QueryBase<T>::QueryBase(const DataSet& dataset)
    : dataset_(dataset)
    , mergeStrategy_(nullptr)
{ }

template<typename T>
inline bool QueryBase<T>::GetNext(T& r)
{
    while (mergeStrategy_->GetNext(r)) {
        if (filterEngine_.Accepts(r))
            return true;
    }
    return false;
}

template<typename T>
inline void QueryIteratorBase<T>::ReadNext(void)
{
    assert(query_);
    if (!query_->GetNext(record_))
        query_ = nullptr;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // QUERYBASE2_H
