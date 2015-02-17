// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#ifndef DICTIONARYBASE_H
#define DICTIONARYBASE_H

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

//
// DictionaryBase is intended to store POD-like structs that are
// uniquely identifiable by one of their fields, referred to as its "key".
//
// Types stored in DictionaryBase must provide read/write access to this
// key field, by providing the following 2 methods:
//
//     // read-only query
//     std::string Key(void) const;
//
//     // write-access to key field & return *this
//     // (T in the signature will be replaced w/ the explicit type)
//     T& Key(const std::string& key);
//
// The implementation of DictionaryBase provides a quick lookup (a la std::map),
// while retaining insertion order for iterating (at the cost of duplicating keys).
//
// Client code must ensure their elements are unique. DictionaryBase will not
// allow elements with identical keys.
//

template<typename T>
class DictionaryBase
{
public:
    typedef std::vector<T>                         DataContainer;
    typedef typename DataContainer::iterator       Iterator;
    typedef typename DataContainer::const_iterator ConstIterator;
    typedef std::map<std::string, size_t>          LookupMap;

// main interface
public:

    // adds value to dictionary, if key is unique
    // returns true if entry was added
    inline bool Add(const T& value);
    inline bool Add(const std::string& key);

    // adds values to dictionary, if keys are unique
    // returns true if all entries added
    inline bool Add(const std::vector<T>& valueList);

    // clear out entries
    inline void Clear(void);

    // returns true if dictionary contains the requested lookup value
    inline bool Contains(const std::string& key) const;
    inline bool Contains(const T& value) const;

    // returns index of key, -1 if not found
    inline int IndexOf(const std::string& key) const;
    inline int IndexOf(const T& value) const;

    // is dictionary empty
    inline bool IsEmpty(void) const;

    // remove requested element
    // returns true if element was actually removed (i.e. returns false if key was not found anyway)
    inline bool Remove(const std::string& key);
    inline bool Remove(const T& value);

    // number of elements
    inline size_t Size(void) const;

    // retrieves a modifiable reference to the object associated with this ID
    // if none exists, a new one is constructed with the associated key and returned
    inline T& operator[](const std::string& key);

    // retrieves a const reference to the object associated with this ID
    // if none exists, throws std::out_of_range (a la std::map::at())
    inline const T& At(const std::string& key) const;

    // retrieves a const reference to the object at index
    // if none exists, throws std::out_of_range
    inline const T& At(const size_t index) const;

// iterators
public:
    inline Iterator      Begin(void);              // returns iterator to begin()
    inline ConstIterator Begin(void) const;        // returns const_iterator to begin()
    inline ConstIterator ConstBegin(void) const;   // returns const_iterator to begin()
    inline Iterator      End(void);                // returns iterator to end()
    inline ConstIterator End(void) const;          // returns const_iterator to end()
    inline ConstIterator ConstEnd(void) const;     // returns const_iterator to end()

// data members
private:
    DataContainer data_;
    LookupMap lookupData_;
};

template<typename T>
inline bool DictionaryBase<T>::Add(const T& value)
{
    if (IsEmpty() || !Contains(value.Key())) {
        data_.push_back(value);
        lookupData_[value.Key()] = data_.size() - 1;
        return true;
    }
    return false;
}

template<typename T>
inline bool DictionaryBase<T>::Add(const std::string& key)
{
    T value;
    value.Key(key);
    return this->Add(value);
}

template<typename T>
inline bool DictionaryBase<T>::Add(const std::vector<T>& valueList)
{
    bool success = true;
    const auto end = valueList.cend();
    for (auto iter = valueList.cbegin(); iter!= end; ++iter)
        success &= this->Add(*iter);
    return success;
}

template<typename T>
inline typename DictionaryBase<T>::Iterator DictionaryBase<T>::Begin(void)
{
    return data_.begin();
}

template<typename T>
inline typename DictionaryBase<T>::ConstIterator DictionaryBase<T>::Begin(void) const
{
    return data_.cbegin();
}

template<typename T>
inline void DictionaryBase<T>::Clear(void)
{
    data_.clear();
    lookupData_.clear();
}

template<typename T>
inline typename DictionaryBase<T>::ConstIterator DictionaryBase<T>::ConstBegin(void) const
{
    return data_.cbegin();
}

template<typename T>
inline typename DictionaryBase<T>::ConstIterator DictionaryBase<T>::ConstEnd(void) const
{
    return data_.cend();
}

template<typename T>
inline bool DictionaryBase<T>::Contains(const std::string& key) const
{
    return (lookupData_.find(key) != lookupData_.cend());
}

template<typename T>
inline bool DictionaryBase<T>::Contains(const T& value) const
{
    return this->Contains(value.Key());
}

template<typename T>
inline typename DictionaryBase<T>::Iterator DictionaryBase<T>::End(void)
{
    return data_.end();
}

template<typename T>
inline typename DictionaryBase<T>::ConstIterator DictionaryBase<T>::End(void) const
{
    return data_.cend();
}

template<typename T>
inline int DictionaryBase<T>::IndexOf(const std::string& key) const
{
    const auto iter = lookupData_.find(key);
    if (iter == lookupData_.cend())
        return -1;
    return iter->second;
}

template<typename T>
inline int DictionaryBase<T>::IndexOf(const T& value) const
{
    return IndexOf(value.Key());
}

template<typename T>
inline bool DictionaryBase<T>::IsEmpty(void) const
{
    return data_.empty();
}

template<typename T>
inline bool DictionaryBase<T>::Remove(const std::string& key)
{
    // skip if empty dictionary or if ID unknown
    if (IsEmpty() || !Contains(key))
        return false;

    // update 'lookup index' for every entry after @readGroupId
    const size_t indexToRemove = lookupData_[key];
    const size_t numEntries = data_.size();
    for (size_t i = indexToRemove+1; i < numEntries; ++i) {
        const T& e = data_.at(i);
        --lookupData_[e.Key()];
    }

    // erase entry from containers
    data_.erase(Begin() + indexToRemove);
    lookupData_.erase(key);
    return true;
}

template<typename T>
inline bool DictionaryBase<T>::Remove(const T& value)
{
    return this->Remove(value.Key());
}

template<typename T>
inline size_t DictionaryBase<T>::Size(void) const
{
    return data_.size();
}

template<typename T>
inline T& DictionaryBase<T>::operator[](const std::string& key)
{
    const auto iter = lookupData_.find(key);
    if (iter == lookupData_.cend()) {
        T newElement;
        newElement.Key(key);
        data_.push_back(newElement);
        lookupData_[key] = data_.size() - 1;
        return data_.at(data_.size()-1);
    }
    return data_.at(iter->second);
}

template<typename T>
inline const T& DictionaryBase<T>::At(const std::string& key) const
{
    const auto iter = lookupData_.find(key);
    if (iter == lookupData_.cend())
        throw std::out_of_range("unknown key");
    return data_.at(iter->second);
}

template<typename T>
inline const T& DictionaryBase<T>::At(const size_t index) const
{
    if ( index >= data_.size() )
        throw std::out_of_range("invalid index");
    return data_.at(index);
}

} // namespace BAM
} // namespace PacBio

#endif // DICTIONARYBASE_H
