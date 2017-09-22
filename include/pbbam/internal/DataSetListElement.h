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

#ifndef DATASETLISTELEMENT_H
#define DATASETLISTELEMENT_H

#include "pbbam/internal/DataSetElement.h"

namespace PacBio {
namespace BAM {
namespace internal {

//
// adds iterators for convenience
//
template<class T> class DataSetListElement;

template<class T>
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

template<class T>
class DataSetListIterator : public DataSetListIteratorBase<T>
{
public:
    DataSetListIterator(const DataSetListElement<T>* parent, size_t i);
    T& operator*();
    T* operator->();
    DataSetListIterator<T>& operator++();
    DataSetListIterator<T> operator++(int);
};

template<class T>
class DataSetListConstIterator : public DataSetListIteratorBase<T>
{
public:
    DataSetListConstIterator(const DataSetListElement<T>* parent, size_t i);
    const T& operator*() const;
    const T* operator->() const;
    DataSetListConstIterator<T>& operator++();
    DataSetListConstIterator<T> operator++(int);
};

template<class T>
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

} // namespace internal
} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/DataSetListElement.inl"

#endif // DATASETLISTELEMENT_H
