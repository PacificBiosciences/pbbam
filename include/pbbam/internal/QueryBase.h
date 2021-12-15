#ifndef PBBAM_QUERYBASE_H
#define PBBAM_QUERYBASE_H

#include <pbbam/Config.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/DataSet.h>

#include <iterator>
#include <memory>
#include <vector>

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

    bool operator==(const QueryIteratorBase<T>& other) const noexcept;
    bool operator!=(const QueryIteratorBase<T>& other) const noexcept;

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
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    QueryIterator() = default;
    QueryIterator(QueryBase<T>& query);

    T& operator*() noexcept;
    T* operator->() noexcept;

    QueryIterator<T>& operator++();
    QueryIterator<T> operator++(int);
};

template <typename T>
class QueryConstIterator : public QueryIteratorBase<T>
{
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = const T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    QueryConstIterator() = default;
    QueryConstIterator(const QueryBase<T>& query);

    const T& operator*() const noexcept;
    const T* operator->() const noexcept;

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

#include <pbbam/internal/QueryBase.inl>

#endif  // PBBAM_QUERYBASE_H
