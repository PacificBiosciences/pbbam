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
// Author: Derek Barnett

#include "pbbam/BamRecord.h"
#include "pbbam/PbiFile.h"
#include "pbbam/PbiIndex.h"
#include "pbbam/PbiRawData.h"

#include <algorithm>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

// --------------------------
// Ordered Lookup Container (e.g. map)
// --------------------------

template<typename T>
class OrderedLookup
{
public:
    typedef T         KeyType;
    typedef IndexList ValueType;
    typedef std::map<KeyType, ValueType> ContainerType;
    typedef typename ContainerType::const_iterator IterType;

public:
    OrderedLookup(void);
    OrderedLookup(size_t n);
    OrderedLookup(const ContainerType& data);
    OrderedLookup(ContainerType&& data);
    OrderedLookup(const std::vector<T>& rawData);
    OrderedLookup(std::vector<T>&& rawData);

public:
    bool operator==(const OrderedLookup<T>& other) const;
    bool operator!=(const OrderedLookup<T>& other) const;

public:
    IndexList LookupIndices(const KeyType& key,
                            const CompareType& compare) const;

private:
    IndexList LookupInclusiveRange(const IterType& begin,
                                   const IterType& end) const;

    IndexList LookupExclusiveRange(const IterType& begin,
                                   const IterType& end,
                                   const KeyType& key) const;

private:
    ContainerType data_;
};

// --------------------------
// Unordered Lookup Container (e.g. hash)
// --------------------------

template<typename T>
class UnorderedLookup
{
public:
    typedef T         KeyType;
    typedef IndexList ValueType;
    typedef std::unordered_map<KeyType, ValueType> ContainerType;

public:
    UnorderedLookup(void);
    UnorderedLookup(size_t n);
    UnorderedLookup(const ContainerType& data);
    UnorderedLookup(ContainerType&& data);
    UnorderedLookup(const std::vector<T>& rawData);
    UnorderedLookup(std::vector<T>&& rawData);

public:
    bool operator==(const UnorderedLookup<T>& other) const;
    bool operator!=(const UnorderedLookup<T>& other) const;

public:
    IndexList LookupIndices(const KeyType& key,
                            const CompareType& compare) const;

private:
    template<typename Compare>
    IndexList LookupHelper(const KeyType& key, const Compare& cmp) const;

private:
    ContainerType data_;
};

// ----------------
// Subread Data
// ----------------

struct BasicLookupData
{
    // ctors
    BasicLookupData(void);
    BasicLookupData(const PbiRawBasicData& rawData);
//    SubreadLookupData(PbiRawSubreadData&& rawData);

    // add offset data to index result blocks
    void ApplyOffsets(IndexResultBlocks& blocks) const;

    template<typename T>
    IndexList Indices(const SubreadField& field,
                      const T& value,
                      const CompareType& compareType) const;

    template<typename T>
    IndexList IndicesMulti(const SubreadField& field,
                           const std::vector<T>& values) const;

    // map ordering doesn't make sense, optimize for direct lookup
    UnorderedLookup<int32_t> rgId_;

    // numeric comparisons make sense, keep key ordering preserved
    OrderedLookup<int32_t> qStart_;
    OrderedLookup<int32_t> qEnd_;
    OrderedLookup<int32_t> holeNumber_;
    OrderedLookup<float>   readQual_;

    // see if this works, or if can use unordered, 'direct' query
    OrderedLookup<uint8_t> ctxtFlag_;

    // offsets
    std::vector<int64_t> fileOffset_;
};

// -----------------
// Mapped Data
// -----------------

struct MappedLookupData
{
    // ctors
    MappedLookupData(void);
    MappedLookupData(const PbiRawMappedData& rawData);
//    MappedLookupData(PbiRawMappedData&& rawData);

    template<typename T>
    IndexList Indices(const MappedField& field,
                      const T& value,
                      const CompareType& compareType) const;

    template<typename T>
    IndexList IndicesMulti(const MappedField& field,
                           const std::vector<T>& values) const;

    // numeric comparisons make sense, keep key ordering preserved
    OrderedLookup<int32_t>  tId_;
    OrderedLookup<uint32_t> tStart_;
    OrderedLookup<uint32_t> tEnd_;
    OrderedLookup<uint32_t> aStart_;
    OrderedLookup<uint32_t> aEnd_;
    OrderedLookup<uint32_t> nM_;
    OrderedLookup<uint32_t> nMM_;
    OrderedLookup<uint8_t>  mapQV_;

    // generated, not stored in PBI
    OrderedLookup<uint32_t> nIns_;
    OrderedLookup<uint32_t> nDel_;

    // no need for map overhead, just store direct indices
    IndexList reverseStrand_;
    IndexList forwardStrand_;
};

// ------------------
// Reference Data
// ------------------

struct ReferenceLookupData
{
    // ctors
    ReferenceLookupData(void);
    ReferenceLookupData(const PbiRawReferenceData& rawData);
//    ReferenceLookupData(PbiRawReferenceData&& rawData);

    IndexRange Indices(const int32_t tId) const;

    // references_[tId] = (begin, end) indices
    // into SubreadLookupData::fileOffset_
    std::unordered_map<int32_t, IndexRange> references_;
};

// ---------------
// Barcode Data
// ---------------

struct BarcodeLookupData
{
    // ctors
    BarcodeLookupData(void);
    BarcodeLookupData(const PbiRawBarcodeData& rawData);
//    BarcodeLookupData(PbiRawBarcodeData&& rawData);

    template<typename T>
    IndexList Indices(const BarcodeField& field,
                      const T& value,
                      const CompareType& compareType) const;

    template<typename T>
    IndexList IndicesMulti(const BarcodeField& field,
                           const std::vector<T>& values) const;

    // numeric comparisons make sense, keep key ordering preserved
    OrderedLookup<uint16_t> bcLeft_;
    OrderedLookup<uint16_t> bcRight_;
    OrderedLookup<uint8_t>  bcQual_;
};

// --------------------------
// Pbi Lookup Aggregate
// --------------------------

class PbiIndexPrivate
{
public:
    PbiIndexPrivate(void);
    PbiIndexPrivate(const PbiRawData& rawIndex);
    PbiIndexPrivate(PbiRawData&& rawIndex);

    std::unique_ptr<PbiIndexPrivate> DeepCopy(void) const;

public:
    bool HasSection(const PbiFile::Section flag) const;
    void SetSection(const PbiFile::Section flag, bool ok = true);

public:

    template<typename T>
    IndexList Indices(const SubreadField& field,
                      const T& value,
                      const CompareType& compareType) const;

    template<typename T>
    IndexList Indices(const MappedField& field,
                      const T& value,
                      const CompareType& compareType) const;

    template<typename T>
    IndexList Indices(const BarcodeField& field,
                      const T& value,
                      const CompareType& compareType) const;

    template<typename T>
    IndexList IndicesMulti(const SubreadField& field,
                           const T& value) const;

    template<typename T>
    IndexList IndicesMulti(const MappedField& field,
                           const T& value) const;

    template<typename T>
    IndexList IndicesMulti(const BarcodeField& field,
                           const T& value) const;

    template<typename T>
    IndexResultBlocks Lookup(const SubreadField& field,
                             const T& value,
                             const CompareType& compareType) const;

    template<typename T>
    IndexResultBlocks Lookup(const MappedField& field,
                             const T& value,
                             const CompareType& compareType) const;

    template<typename T>
    IndexResultBlocks Lookup(const BarcodeField& field,
                             const T& value,
                             const CompareType& compareType) const;

    template<typename T>
    IndexResultBlocks LookupMulti(const SubreadField& field,
                                  const std::vector<T>& values) const;

    template<typename T>
    IndexResultBlocks LookupMulti(const MappedField& field,
                                  const std::vector<T>& values) const;

    template<typename T>
    IndexResultBlocks LookupMulti(const BarcodeField& field,
                                  const std::vector<T>& values) const;

    IndexResultBlocks LookupReference(const int32_t tId) const;

private:
    IndexResultBlocks MergeBlocksWithOffsets(const IndexList& indices) const;

public:
    PbiFile::VersionEnum version_;
    PbiFile::Sections sections_;
    uint32_t numReads_;

    // lookup structures
    BasicLookupData     basicData_;
    MappedLookupData    mappedData_;
    ReferenceLookupData referenceData_;
    BarcodeLookupData   barcodeData_;

private:
    // not-implemented - ensure no copy
    PbiIndexPrivate(const PbiIndexPrivate& other);
    PbiIndexPrivate& operator=(const PbiIndexPrivate& other);
};

// ----------------
// helper methods
// ----------------

inline IndexResultBlocks mergedIndexBlocks(IndexList&& indices)
{
    if (indices.empty())
        return IndexResultBlocks();
    std::sort(indices.begin(), indices.end());

    IndexResultBlocks result;
    result.push_back(IndexResultBlock(indices.at(0), 1));
    const size_t numIndices = indices.size();
    for (size_t i = 1; i < numIndices; ++i) {
        if (indices.at(i) == indices.at(i-1)+1)
            ++result.back().numReads_;
        else
            result.push_back(IndexResultBlock(indices.at(i), 1));
    }
    return result;
}

inline IndexResultBlocks mergedIndexBlocks(const IndexList& indices)
{
    IndexList copy = indices;
    return mergedIndexBlocks(std::move(copy));
}

inline size_t nullIndex(void)
{ return static_cast<size_t>(-1); }

inline
void pushBackIndices(IndexList& result,
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
inline OrderedLookup<T>::OrderedLookup(void) { }

template<typename T>
inline OrderedLookup<T>::OrderedLookup(size_t n)
{ data_.reserve(n); }

template<typename T>
inline OrderedLookup<T>::OrderedLookup(const ContainerType& data)
    : data_(data)
{ }

template<typename T>
inline OrderedLookup<T>::OrderedLookup(ContainerType&& data)
    : data_(std::move(data))
{ }

template<typename T>
inline OrderedLookup<T>::OrderedLookup(const std::vector<T>& rawData)
{
    const size_t numElements = rawData.size();
    for (size_t i = 0; i < numElements; ++i)
        data_[ rawData.at(i) ].push_back(i);
}

template<typename T>
inline OrderedLookup<T>::OrderedLookup(std::vector<T>&& rawData)
{
    const size_t numElements = rawData.size();
    for (size_t i = 0; i < numElements; ++i)
        data_[ rawData.at(i) ].push_back(i);
}

template<typename T>
inline bool OrderedLookup<T>::operator==(const OrderedLookup<T>& other) const
{ return data_ == other.data_; }

template<typename T>
inline bool OrderedLookup<T>::operator!=(const OrderedLookup<T>& other) const
{ return !(*this == other); }

template<typename T>
inline IndexList
OrderedLookup<T>::LookupInclusiveRange(const IterType& begin,
                                       const IterType& end) const
{
    IndexList result;
    for ( auto iter = begin; iter != end; ++iter )
        pushBackIndices(result, iter->second);
    std::sort(result.begin(), result.end());
    return result;
}

template<typename T>
inline IndexList
OrderedLookup<T>::LookupExclusiveRange(const IterType& begin,
                                       const IterType& end,
                                       const KeyType& key) const
{
    IndexList result;
    for ( auto iter = begin; iter != end; ++iter ) {
        if (iter->first != key)
            pushBackIndices(result, iter->second);
    }
    std::sort(result.begin(), result.end());
    return result;
}

template<typename T>
inline IndexList
OrderedLookup<T>::LookupIndices(const OrderedLookup::KeyType& key,
                                const CompareType& compare) const
{
    const IterType begin = data_.cbegin();
    const IterType end   = data_.cend();
    switch(compare)
    {
        case CompareType::EQUAL:
        {
            const auto found = data_.find(key);
            if (found != end)
                return found->second;
            return IndexList();
        }
        case CompareType::LESS_THAN:          return LookupExclusiveRange(begin, data_.upper_bound(key), key);
        case CompareType::LESS_THAN_EQUAL:    return LookupInclusiveRange(begin, data_.upper_bound(key));
        case CompareType::GREATER_THAN:       return LookupExclusiveRange(data_.lower_bound(key), end, key);
        case CompareType::GREATER_THAN_EQUAL: return LookupInclusiveRange(data_.lower_bound(key), end);
        case CompareType::NOT_EQUAL:          return LookupExclusiveRange(begin, end, key);
        default:
            assert(false);
    }
    return IndexList();
}

// -----------------
// UnorderedLookup
// -----------------

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(void) { }

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(size_t n)
{ data_.reserve(n); }

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(const ContainerType& data)
    : data_(data)
{ }

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(ContainerType&& data)
    : data_(std::move(data))
{ }

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(const std::vector<T> &rawData)
{
    const size_t numElements = rawData.size();
    for (size_t i = 0; i < numElements; ++i)
        data_[ rawData.at(i) ].push_back(i);
}

template<typename T>
inline UnorderedLookup<T>::UnorderedLookup(std::vector<T>&& rawData)
{
    const size_t numElements = rawData.size();
    for (size_t i = 0; i < numElements; ++i)
        data_[ rawData.at(i) ].push_back(i);
}

template<typename T>
inline bool UnorderedLookup<T>::operator==(const UnorderedLookup<T>& other) const
{ return data_ == other.data_; }

template<typename T>
inline bool UnorderedLookup<T>::operator!=(const UnorderedLookup<T>& other) const
{ return !(*this == other); }

template<typename T>
template<typename Compare>
inline IndexList
UnorderedLookup<T>::LookupHelper(const UnorderedLookup::KeyType& key,
                                 const Compare& cmp) const
{
    auto iter = data_.cbegin();
    const auto end = data_.cend();
    IndexList result; // init with some avg size ??
    for ( ; iter != end; ++iter ) {
        const auto e = (iter->first);
        if (cmp(e, key))
            pushBackIndices(result, iter->second);
    }
    std::sort(result.begin(), result.end());
    return result;
}

template<typename T>
inline IndexList
UnorderedLookup<T>::LookupIndices(const UnorderedLookup::KeyType& key,
                                  const CompareType& compare) const
{
    switch (compare) {
        case CompareType::EQUAL:
        {
            const auto found = data_.find(key);
            if (found != data_.cend())
                return found->second;
            else
                return IndexList();
        }
        case CompareType::LESS_THAN:          return LookupHelper(key, std::less<KeyType>());
        case CompareType::LESS_THAN_EQUAL:    return LookupHelper(key, std::less_equal<KeyType>());
        case CompareType::GREATER_THAN:       return LookupHelper(key, std::greater<KeyType>());
        case CompareType::GREATER_THAN_EQUAL: return LookupHelper(key, std::greater_equal<KeyType>());
        case CompareType::NOT_EQUAL:          return LookupHelper(key, std::not_equal_to<KeyType>());
        default:
            assert(false);
    }
    return IndexList();
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
inline IndexList BasicLookupData::Indices(const SubreadField& field,
                                            const T& value,
                                            const CompareType& compareType) const
{
    switch(field) {
        case SubreadField::RG_ID:        return rgId_.LookupIndices(value, compareType);
        case SubreadField::Q_START:      return qStart_.LookupIndices(value, compareType);
        case SubreadField::Q_END:        return qEnd_.LookupIndices(value, compareType);
        case SubreadField::ZMW:          return holeNumber_.LookupIndices(value, compareType);
        case SubreadField::READ_QUALITY: return readQual_.LookupIndices(value, compareType);
        case SubreadField::CONTEXT_FLAG: return ctxtFlag_.LookupIndices(value, compareType);

        case SubreadField::VIRTUAL_OFFSET : // fall-through, not supported this way
        default:
            assert(false);
    }
    return IndexList();
}

template<typename T>
inline IndexList BasicLookupData::IndicesMulti(const SubreadField& field,
                                                 const std::vector<T>& values) const
{
    IndexList result;
    for (auto value : values) {
        const IndexList& valueIndices = Indices(field, value, CompareType::EQUAL);
        result.reserve(result.size() + valueIndices.size());
        for (auto i : valueIndices)
            result.push_back(i);
    }
    return result;
}

// -------------------
// MappedLookupData
// -------------------

template<typename T>
inline IndexList MappedLookupData::Indices(const MappedField& field,
                                           const T& value,
                                           const CompareType& compareType) const
{
    switch(field) {
        case MappedField::T_ID:        return tId_.LookupIndices(value, compareType);
        case MappedField::T_START:     return tStart_.LookupIndices(value, compareType);
        case MappedField::T_END:       return tEnd_.LookupIndices(value, compareType);
        case MappedField::A_START:     return aStart_.LookupIndices(value, compareType);
        case MappedField::A_END:       return aEnd_.LookupIndices(value, compareType);
        case MappedField::N_M:         return nM_.LookupIndices(value, compareType);
        case MappedField::N_MM:        return nM_.LookupIndices(value, compareType);
        case MappedField::MAP_QUALITY: return mapQV_.LookupIndices(value, compareType);

        // MappedField::STRAND has its own specialization

        default:
            assert(false);
    }
    return IndexList();
}

template<>
inline IndexList MappedLookupData::Indices(const MappedField& field,
                                           const Strand& strand,
                                           const CompareType& compareType) const
{
    assert(field == MappedField::STRAND);

    if (compareType == CompareType::EQUAL) {
        if (strand == Strand::FORWARD)
            return forwardStrand_;
        else
            return reverseStrand_;
    } else if (compareType == CompareType::NOT_EQUAL) {
        if (strand == Strand::FORWARD)
            return reverseStrand_;
        else
            return forwardStrand_;
    }

    // only EQUAL/NOT_EQUAL supported
    assert(false);
    return IndexList();
}

template<typename T>
inline IndexList MappedLookupData::IndicesMulti(const MappedField& field,
                                                const std::vector<T>& values) const
{
    IndexList result;
    for (auto value : values) {
        const IndexList& valueIndices = Indices(field, value, CompareType::EQUAL);
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
        return IndexRange(nullIndex(), nullIndex());
    return found->second;
}

// -------------------
// BarcodeLookupData
// -------------------

template<typename T>
inline IndexList BarcodeLookupData::Indices(const BarcodeField& field,
                                            const T& value,
                                            const CompareType& compareType) const
{
    switch(field) {
        case BarcodeField::BC_LEFT:      return bcLeft_.LookupIndices(value, compareType);
        case BarcodeField::BC_RIGHT:     return bcRight_.LookupIndices(value, compareType);
        case BarcodeField::BC_QUALITY:   return bcQual_.LookupIndices(value, compareType);
        default:
            assert(false);
    }
    return IndexList();
}

template<typename T>
inline IndexList BarcodeLookupData::IndicesMulti(const BarcodeField& field,
                                                 const std::vector<T>& values) const
{
    IndexList result;
    for (auto value : values) {
        const IndexList& valueIndices = Indices(field, value, CompareType::EQUAL);
        result.reserve(result.size() + valueIndices.size());
        for (auto i : valueIndices)
            result.push_back(i);
    }
    return result;
}


// -----------------
// PbiIndexPrivate
// -----------------

inline bool PbiIndexPrivate::HasSection(const PbiFile::Section flag) const
{ return (sections_ & flag) != 0; }

inline void PbiIndexPrivate::SetSection(const PbiFile::Section flag, bool ok)
{ if (ok) sections_ |= flag; else sections_ &= ~flag; }

template<typename T>
inline IndexList
PbiIndexPrivate::Indices(const SubreadField& field,
                         const T& value,
                         const CompareType& compareType) const
{ return basicData_.Indices(field, value, compareType); }

template<typename T>
inline IndexList
PbiIndexPrivate::Indices(const MappedField& field,
                         const T& value,
                         const CompareType& compareType) const
{ return mappedData_.Indices(field, value, compareType); }

template<typename T>
inline IndexList
PbiIndexPrivate::Indices(const BarcodeField& field,
                         const T& value,
                         const CompareType& compareType) const
{ return barcodeData_.Indices(field, value, compareType); }

template<typename T>
inline IndexList
PbiIndexPrivate::IndicesMulti(const SubreadField& field,
                              const T& value) const
{ return basicData_.IndicesMulti(field, value); }

template<typename T>
inline IndexList
PbiIndexPrivate::IndicesMulti(const MappedField& field,
                              const T& value) const
{ return mappedData_.IndicesMulti(field, value); }

template<typename T>
inline IndexList
PbiIndexPrivate::IndicesMulti(const BarcodeField& field,
                              const T& value) const
{ return barcodeData_.IndicesMulti(field, value); }

template<typename T>
inline IndexResultBlocks
PbiIndexPrivate::Lookup(const SubreadField& field,
                        const T& value,
                        const CompareType& compareType) const
{ return MergeBlocksWithOffsets(basicData_.Indices(field, value, compareType)); }

template<typename T>
inline IndexResultBlocks
PbiIndexPrivate::Lookup(const MappedField& field,
                        const T& value,
                        const CompareType& compareType) const
{
    if (!HasSection(PbiFile::MAPPED))
        return IndexResultBlocks();
    return MergeBlocksWithOffsets(mappedData_.Indices(field, value, compareType));
}

template<typename T>
inline IndexResultBlocks
PbiIndexPrivate::Lookup(const BarcodeField& field,
                        const T& value,
                        const CompareType& compareType) const
{
    if (!HasSection(PbiFile::BARCODE))
        return IndexResultBlocks();
    return MergeBlocksWithOffsets(barcodeData_.Indices(field, value, compareType));
}

template<typename T>
inline IndexResultBlocks
PbiIndexPrivate::LookupMulti(const SubreadField& field,
                             const std::vector<T>& values) const
{ return MergeBlocksWithOffsets(basicData_.IndicesMulti(field, values)); }

template<typename T>
inline IndexResultBlocks
PbiIndexPrivate::LookupMulti(const MappedField& field,
                             const std::vector<T>& values) const
{ return MergeBlocksWithOffsets(mappedData_.IndicesMulti(field, values)); }

template<typename T>
inline IndexResultBlocks
PbiIndexPrivate::LookupMulti(const BarcodeField& field,
                             const std::vector<T>& values) const
{ return MergeBlocksWithOffsets(barcodeData_.IndicesMulti(field, values)); }

inline IndexResultBlocks
PbiIndexPrivate::LookupReference(const int32_t tId) const
{
    if (!HasSection(PbiFile::REFERENCE))
        return IndexResultBlocks();
    const IndexRange& indexRange = referenceData_.Indices(tId);
    if (indexRange.first == nullIndex() && indexRange.second == nullIndex())
        return IndexResultBlocks();
    const size_t numReads = indexRange.second - indexRange.first;
    IndexResultBlocks blocks(1, IndexResultBlock(indexRange.first, numReads));
    basicData_.ApplyOffsets(blocks);
    return blocks;
}

inline IndexResultBlocks
PbiIndexPrivate::MergeBlocksWithOffsets(const IndexList& indices) const
{
    IndexResultBlocks blocks = mergedIndexBlocks(indices);
    basicData_.ApplyOffsets(blocks);
    return blocks;
}

} // namespace internal

template<typename FieldType, typename ValueType>
inline IndexRequestBase<FieldType, ValueType>::IndexRequestBase(const FieldType field,
                                                                const ValueType& value,
                                                                const CompareType compareType)
    : field_(field)
    , value_(value)
    , compareType_(compareType)
{ }

template<typename FieldType, typename ValueType>
inline IndexMultiRequestBase<FieldType, ValueType>::IndexMultiRequestBase(const FieldType field,
                                                                          const std::vector<ValueType>& values)
    : field_(field)
    , values_(values)
{ }

template<SubreadField field, typename ValueType>
inline SubreadIndexRequest<field, ValueType>::SubreadIndexRequest(const ValueType& value,
                                                                  const CompareType& compareType)
    : IndexRequestBase<SubreadField, ValueType>(field, value, compareType)
{ }

template<SubreadField field, typename ValueType>
inline SubreadIndexMultiRequest<field, ValueType>::SubreadIndexMultiRequest(const std::vector<ValueType>& values)
    : IndexMultiRequestBase<SubreadField, ValueType>(field, values)
{ }

template<MappedField field, typename ValueType>
inline MappedIndexRequest<field, ValueType>::MappedIndexRequest(const ValueType& value,
                                                                const CompareType& compareType)
    : IndexRequestBase<MappedField, ValueType>(field, value, compareType)
{ }

template<MappedField field, typename ValueType>
inline MappedIndexMultiRequest<field, ValueType>::MappedIndexMultiRequest(const std::vector<ValueType>& values)
    : IndexMultiRequestBase<MappedField, ValueType>(field, values)
{ }

template<BarcodeField field, typename ValueType>
inline BarcodeIndexRequest<field, ValueType>::BarcodeIndexRequest(const ValueType& value,
                                                                  const CompareType& compareType)
    : IndexRequestBase<BarcodeField, ValueType>(field, value, compareType)
{ }

template<BarcodeField field, typename ValueType>
inline BarcodeIndexMultiRequest<field, ValueType>::BarcodeIndexMultiRequest(const std::vector<ValueType>& values)
    : IndexMultiRequestBase<BarcodeField, ValueType>(field, values)
{ }

template <typename FieldType, typename ValueType>
inline IndexList
PbiIndex::RawIndices(const IndexRequestBase<FieldType, ValueType>& request) const
{ return d_->Indices(request.field_, request.value_, request.compareType_); }

template <typename FieldType, typename ValueType>
inline IndexList
PbiIndex::RawIndices(const IndexMultiRequestBase<FieldType, ValueType>& request) const
{ return d_->Indices(request.field_, request.values_); }

template <typename FieldType, typename ValueType>
inline IndexResultBlocks
PbiIndex::Lookup(const IndexRequestBase<FieldType, ValueType>& request) const
{ return d_->Lookup(request.field_, request.value_, request.compareType_); }

template <typename FieldType, typename ValueType>
inline IndexResultBlocks
PbiIndex::Lookup(const IndexMultiRequestBase<FieldType, ValueType>& request) const
{ return d_->LookupMulti(request.field_, request.values_); }

inline IndexResultBlocks PbiIndex::LookupReference(const int32_t tId) const
{ return d_->LookupReference(tId); }

} // namespace BAM
} // namespace PacBio

