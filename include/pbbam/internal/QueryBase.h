// Author: Derek Barnett

#ifndef QUERYBASE_H
#define QUERYBASE_H

#include <cassert>
#include <memory>
#include <vector>
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T>
class QueryBase;

template <typename T>
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

template <typename T>
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

template <typename T>
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

template <typename T>
class QueryBase
{

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
    virtual bool GetNext(T& r) = 0;

protected:
    QueryBase() = default;
};

using IQuery = QueryBase<BamRecord>;
using IGroupQuery = QueryBase<std::vector<BamRecord>>;

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/QueryBase.inl"

#endif  // QUERYBASE_H
