#ifndef PBBAM_QUERYBASE_INL
#define PBBAM_QUERYBASE_INL

#include <pbbam/internal/QueryBase.h>

#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

// -------------------
// QueryIteratorBase
// -------------------

template <typename T>
QueryIteratorBase<T>::QueryIteratorBase(QueryBase<T>& query) : query_{&query}
{
    ReadNext();
}

template <typename T>
bool QueryIteratorBase<T>::operator==(const QueryIteratorBase<T>& other) const noexcept
{
    return query_ == other.query_;
}

template <typename T>
bool QueryIteratorBase<T>::operator!=(const QueryIteratorBase<T>& other) const noexcept
{
    return !(*this == other);
}

// -------------------
// QueryIterator
// -------------------

template <typename T>
QueryIterator<T>::QueryIterator(QueryBase<T>& query) : QueryIteratorBase<T>{query}
{
}

template <typename T>
T& QueryIterator<T>::operator*() noexcept
{
    return QueryIteratorBase<T>::record_;
}

template <typename T>
T* QueryIterator<T>::operator->() noexcept
{
    return &(operator*());
}

template <typename T>
QueryIterator<T>& QueryIterator<T>::operator++()
{
    QueryIteratorBase<T>::ReadNext();
    return *this;
}

template <typename T>
QueryIterator<T> QueryIterator<T>::operator++(int)
{
    QueryIterator<T> result(*this);
    ++(*this);
    return result;
}

// --------------------
// QueryConstIterator
// --------------------

template <typename T>
QueryConstIterator<T>::QueryConstIterator(const QueryBase<T>& query)
    : QueryIteratorBase<T>{const_cast<QueryBase<T>&>(query)}
{
}

template <typename T>
const T& QueryConstIterator<T>::operator*() const noexcept
{
    return QueryIteratorBase<T>::record_;
}

template <typename T>
const T* QueryConstIterator<T>::operator->() const noexcept
{
    return &(operator*());
}

template <typename T>
QueryConstIterator<T>& QueryConstIterator<T>::operator++()
{
    QueryIteratorBase<T>::ReadNext();
    return *this;
}

template <typename T>
QueryConstIterator<T> QueryConstIterator<T>::operator++(int)
{
    QueryConstIterator<T> result(*this);
    ++(*this);
    return result;
}

// -----------
// QueryBase
// -----------

template <typename T>
QueryConstIterator<T> QueryBase<T>::begin() const
{
    return QueryConstIterator<T>(*this);
}

template <typename T>
QueryIterator<T> QueryBase<T>::begin()
{
    return QueryIterator<T>(*this);
}

template <typename T>
QueryConstIterator<T> QueryBase<T>::cbegin() const
{
    return QueryConstIterator<T>(*this);
}

template <typename T>
QueryConstIterator<T> QueryBase<T>::cend() const
{
    return QueryConstIterator<T>();
}

template <typename T>
QueryConstIterator<T> QueryBase<T>::end() const
{
    return QueryConstIterator<T>();
}

template <typename T>
QueryIterator<T> QueryBase<T>::end()
{
    return QueryIterator<T>();
}

template <typename T>
void QueryIteratorBase<T>::ReadNext()
{
    assert(query_);
    if (!query_->GetNext(record_)) {
        query_ = nullptr;
    }
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif // PBBAM_QUERYBASE_INL
