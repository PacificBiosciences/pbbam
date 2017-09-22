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
//
// File Description
/// \file PbiLookupData.inl
/// \brief Inline implementations for the classes used for PBI data lookup.
//
// Author: Derek Barnett

#include "pbbam/PbiLookupData.h"
#include "pbbam/PbiRawData.h"
#include "pbbam/Strand.h"
#include <algorithm>
#include <unordered_set>
#include <cassert>

namespace PacBio {
namespace BAM {

// ----------------
// helper methods
// ----------------

inline IndexResultBlocks mergedIndexBlocks(IndexList&& indices)
{
    if (indices.empty())
        return IndexResultBlocks{ };

    std::sort(indices.begin(), indices.end());
    auto newEndIter = std::unique(indices.begin(), indices.end());
    auto numIndices = std::distance(indices.begin(), newEndIter);
    assert(!indices.empty());
    auto result = IndexResultBlocks{ IndexResultBlock(indices.at(0), 1) };
    for (auto i = 1; i < numIndices; ++i) {
        if (indices.at(i) == indices.at(i-1)+1)
            ++result.back().numReads_;
        else
            result.emplace_back(indices.at(i), 1);
    }
    return result;
}

inline IndexResultBlocks mergedIndexBlocks(const IndexList& indices)
{
    auto copy = indices;
    return mergedIndexBlocks(std::move(copy));
}

inline size_t nullIndex()
{ return static_cast<size_t>(-1); }

inline void pushBackIndices(IndexList& result,
                            const IndexList& toAppend)
{
    result.reserve(result.size() + toAppend.size());
    for (auto element : toAppend)
        result.push_back(element);
}

// -----------------
// OrderedLookup
// -----------------

template<typename T>
inline OrderedLookup<T>::OrderedLookup(const container_type& data)
    : data_(data)
{ }

template<typename T>
inline OrderedLookup<T>::OrderedLookup(container_type&& data)
    : data_(std::move(data))
{ }

template<typename T>
inline OrderedLookup<T>::OrderedLookup(const std::vector<T>& rawData)
{
    const auto numElements = rawData.size();
    for (auto i = decltype(numElements){0}; i < numElements; ++i)
        data_[rawData.at(i)].push_back(i);
}

template<typename T>
inline OrderedLookup<T>::OrderedLookup(std::vector<T>&& rawData)
{
    const auto numElements = rawData.size();
    for (auto i = decltype(numElements){0}; i < numElements; ++i)
        data_[rawData.at(i)].push_back(i);
}

template<typename T>
inline bool OrderedLookup<T>::operator==(const OrderedLookup<T>& other) const
{ return data_ == other.data_; }

template<typename T>
inline bool OrderedLookup<T>::operator!=(const OrderedLookup<T>& other) const
{ return !(*this == other); }

template<typename T>
inline typename OrderedLookup<T>::iterator OrderedLookup<T>::begin()
{ return data_.begin(); }

template<typename T>
inline typename OrderedLookup<T>::const_iterator OrderedLookup<T>::begin() const
{ return data_.cbegin(); }

template<typename T>
inline typename OrderedLookup<T>::const_iterator OrderedLookup<T>::cbegin() const
{ return data_.cbegin(); }

template<typename T>
inline typename OrderedLookup<T>::iterator OrderedLookup<T>::end()
{ return data_.end(); }

template<typename T>
inline typename OrderedLookup<T>::const_iterator OrderedLookup<T>::end() const
{ return data_.cend(); }

template<typename T>
inline typename OrderedLookup<T>::const_iterator OrderedLookup<T>::cend() const
{ return data_.cend(); }

template<typename T>
inline bool OrderedLookup<T>::empty() const
{ return data_.empty(); }

template<typename T>
inline size_t OrderedLookup<T>::size() const
{ return data_.size(); }

template<typename T>
inline IndexList
OrderedLookup<T>::LookupInclusiveRange(const const_iterator &begin,
                                       const const_iterator &end) const
{
    auto result = IndexList{ };
    for (auto iter = begin; iter != end; ++iter)
        pushBackIndices(result, iter->second);
    std::sort(result.begin(), result.end());
    return result;
}

template<typename T>
inline IndexList
OrderedLookup<T>::LookupExclusiveRange(const const_iterator& begin,
                                       const const_iterator& end,
                                       const key_type& key) const
{
    auto result = IndexList{ };
    for (auto iter = begin; iter != end; ++iter) {
        if (iter->first != key)
            pushBackIndices(result, iter->second);
    }
    std::sort(result.begin(), result.end());
    return result;
}

template<typename T>
inline IndexList
OrderedLookup<T>::LookupIndices(const OrderedLookup::key_type& key,
                                const Compare::Type& compare) const
{
    auto begin = data_.cbegin();
    auto end   = data_.cend();
    switch(compare)
    {
        case Compare::EQUAL:
        {
            const auto found = data_.find(key);
            if (found != end)
                return found->second;
            return IndexList();
        }
        case Compare::LESS_THAN:          return LookupExclusiveRange(begin, data_.upper_bound(key), key);
        case Compare::LESS_THAN_EQUAL:    return LookupInclusiveRange(begin, data_.upper_bound(key));
        case Compare::GREATER_THAN:       return LookupExclusiveRange(data_.lower_bound(key), end, key);
        case Compare::GREATER_THAN_EQUAL: return LookupInclusiveRange(data_.lower_bound(key), end);
        case Compare::NOT_EQUAL:          return LookupExclusiveRange(begin, end, key);
        default:
            assert(false);
    }
    return IndexList{ };
}

template<typename T>
inline std::vector<T> OrderedLookup<T>::Unpack() const
{
    auto result = std::vector<T>{ };
    auto iter = cbegin();
    const auto end = cend();
    for ( ; iter != end; ++iter ) {
        const auto& indices = iter->second;
        for (auto&& i : indices) {
            if (result.size() <= i)
                result.resize(i+1);
            result[i] = iter->first;
        }
    }
    return result;
}

// -----------------
// UnorderedLookup
// -----------------

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(const container_type& data)
    : data_(data)
{ }

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(container_type&& data)
    : data_(std::move(data))
{ }

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(const std::vector<T>& rawData)
{
    const auto numElements = rawData.size();
    for (auto i = decltype(numElements){0}; i < numElements; ++i)
        data_[rawData.at(i)].push_back(i);
}

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(std::vector<T>&& rawData)
{
    const auto numElements = rawData.size();
    for (auto i = decltype(numElements){0}; i < numElements; ++i)
        data_[rawData.at(i)].push_back(i);
}

template<typename T>
inline bool UnorderedLookup<T>::operator==(const UnorderedLookup<T>& other) const
{ return data_ == other.data_; }

template<typename T>
inline bool UnorderedLookup<T>::operator!=(const UnorderedLookup<T>& other) const
{ return !(*this == other); }

template<typename T>
inline typename UnorderedLookup<T>::iterator UnorderedLookup<T>::begin()
{ return data_.begin(); }

template<typename T>
inline typename UnorderedLookup<T>::const_iterator UnorderedLookup<T>::begin() const
{ return data_.cbegin(); }

template<typename T>
inline typename UnorderedLookup<T>::const_iterator UnorderedLookup<T>::cbegin() const
{ return data_.cbegin(); }

template<typename T>
inline typename UnorderedLookup<T>::iterator UnorderedLookup<T>::end()
{ return data_.end(); }

template<typename T>
inline typename UnorderedLookup<T>::const_iterator UnorderedLookup<T>::end() const
{ return data_.cend(); }

template<typename T>
inline typename UnorderedLookup<T>::const_iterator UnorderedLookup<T>::cend() const
{ return data_.cend(); }

template<typename T>
inline bool UnorderedLookup<T>::empty() const
{ return data_.empty(); }

template<typename T>
inline size_t UnorderedLookup<T>::size() const
{ return data_.size(); }

template<typename T>
template<typename Compare>
inline IndexList
UnorderedLookup<T>::LookupHelper(const UnorderedLookup::key_type& key,
                                 const Compare& cmp) const
{
    auto result = IndexList{ }; // init with some avg size ??
    const auto end = data_.cend();
    for (auto iter = data_.cbegin(); iter != end; ++iter) {
        const auto e = (iter->first);
        if (cmp(e, key))
            pushBackIndices(result, iter->second);
    }
    std::sort(result.begin(), result.end());
    return result;
}

template<typename T>
inline IndexList
UnorderedLookup<T>::LookupIndices(const UnorderedLookup::key_type& key,
                                  const Compare::Type& compare) const
{
    switch (compare) {
        case Compare::EQUAL:
        {
            const auto found = data_.find(key);
            if (found != data_.cend())
                return found->second;
            else
                return IndexList();
        }
        case Compare::LESS_THAN:          return LookupHelper(key, std::less<key_type>());
        case Compare::LESS_THAN_EQUAL:    return LookupHelper(key, std::less_equal<key_type>());
        case Compare::GREATER_THAN:       return LookupHelper(key, std::greater<key_type>());
        case Compare::GREATER_THAN_EQUAL: return LookupHelper(key, std::greater_equal<key_type>());
        case Compare::NOT_EQUAL:          return LookupHelper(key, std::not_equal_to<key_type>());
        default:
            assert(false);
    }
    return IndexList{ };
}

template<typename T>
inline std::vector<T> UnorderedLookup<T>::Unpack() const
{
    auto result = std::vector<T>{ };
    auto iter = cbegin();
    const auto end = cend();
    for ( ; iter != end; ++iter ) {
        const auto& indices = iter->second;
        for (auto&& i : indices) {
            if (result.size() <= i)
                result.resize(i+1);
            result[i] = iter->first;
        }
    }
    return result;
}

// -------------------
// SubreadLookupData
// -------------------

inline
void BasicLookupData::ApplyOffsets(IndexResultBlocks& blocks) const
{
    for (IndexResultBlock& block : blocks)
        block.virtualOffset_ = fileOffset_.at(block.firstIndex_);
}

template<typename T>
inline IndexList BasicLookupData::Indices(const BasicLookupData::Field& field,
                                            const T& value,
                                            const Compare::Type& compareType) const
{
    switch(field) {
        case BasicLookupData::RG_ID:        return rgId_.LookupIndices(value, compareType);
        case BasicLookupData::Q_START:      return qStart_.LookupIndices(value, compareType);
        case BasicLookupData::Q_END:        return qEnd_.LookupIndices(value, compareType);
        case BasicLookupData::ZMW:          return holeNumber_.LookupIndices(value, compareType);
        case BasicLookupData::READ_QUALITY: return readQual_.LookupIndices(value, compareType);
        case BasicLookupData::CONTEXT_FLAG: return ctxtFlag_.LookupIndices(value, compareType);

        case BasicLookupData::VIRTUAL_OFFSET : // fall-through, not supported this way
        default:
            assert(false);
    }
    return IndexList{ };
}

template<typename T>
inline IndexList BasicLookupData::IndicesMulti(const BasicLookupData::Field& field,
                                                 const std::vector<T>& values) const
{
    auto result = IndexList{ };
    for (auto value : values) {
        const auto valueIndices = Indices(field, value, Compare::EQUAL);
        result.reserve(result.size() + valueIndices.size());
        for (auto i : valueIndices)
            result.push_back(i);
    }
    return result;
}

inline const std::vector<int64_t>& BasicLookupData::VirtualFileOffsets() const
{ return fileOffset_; }

// -------------------
// MappedLookupData
// -------------------

template<typename T>
inline IndexList MappedLookupData::Indices(const MappedLookupData::Field& field,
                                           const T& value,
                                           const Compare::Type& compareType) const
{
    switch(field) {
        case MappedLookupData::T_ID:        return tId_.LookupIndices(value, compareType);
        case MappedLookupData::T_START:     return tStart_.LookupIndices(value, compareType);
        case MappedLookupData::T_END:       return tEnd_.LookupIndices(value, compareType);
        case MappedLookupData::A_START:     return aStart_.LookupIndices(value, compareType);
        case MappedLookupData::A_END:       return aEnd_.LookupIndices(value, compareType);
        case MappedLookupData::N_M:         return nM_.LookupIndices(value, compareType);
        case MappedLookupData::N_MM:        return nMM_.LookupIndices(value, compareType);
        case MappedLookupData::N_DEL:       return nDel_.LookupIndices(value, compareType);
        case MappedLookupData::N_INS:       return nIns_.LookupIndices(value, compareType);
        case MappedLookupData::MAP_QUALITY: return mapQV_.LookupIndices(value, compareType);

        // MappedField::STRAND has its own specialization

        default:
            assert(false);
    }
    return IndexList{ };
}

template<>
inline IndexList MappedLookupData::Indices(const MappedLookupData::Field& field,
                                           const Strand& strand,
                                           const Compare::Type& compareType) const
{
    assert(field == MappedLookupData::STRAND);
//    ()field; // quash warnings building in release mode

    if (compareType == Compare::EQUAL) {
        if (strand == Strand::FORWARD)
            return forwardStrand_;
        else
            return reverseStrand_;
    } else if (compareType == Compare::NOT_EQUAL) {
        if (strand == Strand::FORWARD)
            return reverseStrand_;
        else
            return forwardStrand_;
    }

    // only EQUAL/NOT_EQUAL supported
    assert(false);
    return IndexList{ };
}

template<typename T>
inline IndexList MappedLookupData::IndicesMulti(const MappedLookupData::Field& field,
                                                const std::vector<T>& values) const
{
    auto result = IndexList{ };
    for (auto value : values) {
        auto valueIndices = Indices(field, value, Compare::EQUAL);
        result.reserve(result.size() + valueIndices.size());
        for (auto i : valueIndices)
            result.push_back(i);
    }
    return result;
}


// ---------------------
// ReferenceLookupData
// ---------------------

inline IndexRange ReferenceLookupData::Indices(const int32_t tId) const
{
    auto found = references_.find(tId);
    if (found == references_.cend())
        return IndexRange{ nullIndex(), nullIndex() };
    return found->second;
}

// -------------------
// BarcodeLookupData
// -------------------

template<typename T>
inline IndexList BarcodeLookupData::Indices(const BarcodeLookupData::Field &field,
                                            const T& value,
                                            const Compare::Type &compareType) const
{
    switch(field) {
        case BarcodeLookupData::BC_FORWARD:      return bcForward_.LookupIndices(value, compareType);
        case BarcodeLookupData::BC_REVERSE:     return bcReverse_.LookupIndices(value, compareType);
        case BarcodeLookupData::BC_QUALITY:   return bcQual_.LookupIndices(value, compareType);
        default:
            assert(false);
    }
    return IndexList{ };
}

template<typename T>
inline IndexList BarcodeLookupData::IndicesMulti(const BarcodeLookupData::Field &field,
                                                 const std::vector<T>& values) const
{
    IndexList result;
    for (auto value : values) {
        const IndexList& valueIndices = Indices(field, value, Compare::EQUAL);
        result.reserve(result.size() + valueIndices.size());
        for (auto i : valueIndices)
            result.push_back(i);
    }
    return result;
}

} // namespace BAM
} // namespace PacBio
