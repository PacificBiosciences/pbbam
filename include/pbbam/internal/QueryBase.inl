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
    : query_{&query}
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
    : QueryIteratorBase<T>{query}
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
    : QueryIteratorBase<T>{const_cast<QueryBase<T>&>(query)}
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
