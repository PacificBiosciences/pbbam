// Author: Derek Barnett

#ifndef DATASETLISTELEMENT_H
#define DATASETLISTELEMENT_H

#include "pbbam/internal/DataSetElement.h"

namespace PacBio {
namespace BAM {
namespace internal {

//
// adds iterators for convenience
//
template <class T>
class DataSetListElement;

template <class T>
class DataSetListIteratorBase
{
public:
    bool operator==(const DataSetListIteratorBase<T>& other) const;
    bool operator!=(const DataSetListIteratorBase<T>& other) const;

protected:
    DataSetListIteratorBase(const DataSetListElement<T>* parent, size_t i);
    void ReadNext();

protected:
    const DataSetListElement<T>* parent_;
    size_t index_;
};

template <class T>
class DataSetListIterator : public DataSetListIteratorBase<T>
{
public:
    DataSetListIterator(const DataSetListElement<T>* parent, size_t i);
    T& operator*();
    T* operator->();
    DataSetListIterator<T>& operator++();
    DataSetListIterator<T> operator++(int);
};

template <class T>
class DataSetListConstIterator : public DataSetListIteratorBase<T>
{
public:
    DataSetListConstIterator(const DataSetListElement<T>* parent, size_t i);
    const T& operator*() const;
    const T* operator->() const;
    DataSetListConstIterator<T>& operator++();
    DataSetListConstIterator<T> operator++(int);
};

template <class T>
class DataSetListElement : public DataSetElement
{
public:
    DataSetListElement(const std::string& label, const XsdType& xsd = XsdType::NONE);

    // child access through index
public:
    const T& operator[](size_t index) const;
    T& operator[](size_t index);
    size_t Size() const;

    // child access through iterators
public:
    DataSetListIterator<T> begin();
    DataSetListConstIterator<T> begin() const;
    DataSetListConstIterator<T> cbegin() const;
    DataSetListIterator<T> end();
    DataSetListConstIterator<T> end() const;
    DataSetListConstIterator<T> cend() const;
};

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/DataSetListElement.inl"

#endif  // DATASETLISTELEMENT_H
