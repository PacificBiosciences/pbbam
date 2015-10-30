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

#ifndef PBILOOKUPDATA_H
#define PBILOOKUPDATA_H

#include "pbbam/Config.h"
#include "pbbam/Compare.h"
#include "pbbam/PbiBasicTypes.h"
#include <deque>
#include <map>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace BAM {

class PbiRawBarcodeData;
class PbiRawBasicData;
class PbiRawMappedData;
class PbiRawReferenceData;

// --------------------------
// Ordered Lookup Container (e.g. map)
// --------------------------

// stores (value, indexlist) pairs, ordered by value (likely using std::map)
template<typename T>
class OrderedLookup
{
public:
    typedef T                                       key_type;
    typedef IndexList                               value_type;
    typedef std::map<key_type, value_type>          container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

public:
    OrderedLookup(void);
    OrderedLookup(const container_type& data);
    OrderedLookup(container_type&& data);
    OrderedLookup(const std::vector<T>& rawData);
    OrderedLookup(std::vector<T>&& rawData);

public:
    bool operator==(const OrderedLookup<T>& other) const;
    bool operator!=(const OrderedLookup<T>& other) const;

public:
    // STL-compatibility
    iterator begin(void);
    const_iterator begin(void) const;
    const_iterator cbegin(void) const;

    iterator end(void);
    const_iterator end(void) const;
    const_iterator cend(void) const;

    bool empty(void) const;
    size_t size(void) const;

public:
    IndexList LookupIndices(const key_type& key,
                            const Compare::Type& compare) const;

    std::vector<T> Unpack(void) const;

private:
    IndexList LookupInclusiveRange(const const_iterator& begin,
                                   const const_iterator& end) const;

    IndexList LookupExclusiveRange(const const_iterator& begin,
                                   const const_iterator& end,
                                   const key_type& key) const;

private:
    container_type data_;
};

// --------------------------
// Unordered Lookup Container (e.g. hash)
// --------------------------

// stores (value, indexlist) pairs, not ordered (likely using std::unordered_map)
template<typename T>
class UnorderedLookup
{
public:
    typedef T                                        key_type;
    typedef IndexList                                value_type;
    typedef std::unordered_map<key_type, value_type> container_type;
    typedef typename container_type::iterator        iterator;
    typedef typename container_type::const_iterator  const_iterator;

public:
    UnorderedLookup(void);
    UnorderedLookup(const container_type& data);
    UnorderedLookup(container_type&& data);
    UnorderedLookup(const std::vector<T>& rawData);
    UnorderedLookup(std::vector<T>&& rawData);

public:
    bool operator==(const UnorderedLookup<T>& other) const;
    bool operator!=(const UnorderedLookup<T>& other) const;

public:
    iterator begin(void);
    const_iterator begin(void) const;
    const_iterator cbegin(void) const;

    iterator end(void);
    const_iterator end(void) const;
    const_iterator cend(void) const;

    bool empty(void) const;
    size_t size(void) const;

public:
    IndexList LookupIndices(const key_type& key,
                            const Compare::Type& compare) const;

    std::vector<T> Unpack(void) const;

private:
    template<typename Compare>
    IndexList LookupHelper(const key_type& key, const Compare& cmp) const;

private:
    container_type data_;
};

// ----------------
// Subread Data
// ----------------

class PBBAM_EXPORT BasicLookupData
{
public:
    enum Field
    {
        RG_ID
      , Q_START
      , Q_END
      , ZMW
      , READ_QUALITY
      , CONTEXT_FLAG
      , VIRTUAL_OFFSET
    };

public:
    BasicLookupData(void);
    BasicLookupData(const PbiRawBasicData& rawData);

public:

    /// For each block in \p blocks, add its virtual offset.
    /// Readers will use this to seek to block.
    ///
    void ApplyOffsets(IndexResultBlocks& blocks) const;

    /// Helper method that dispatches a filter request to the proper data member,
    /// respecting compareType.
    ///
    /// Client code (e.g. PBI filters) should use this when possible, only touching
    /// raw fields for more complex operations.
    ///
    template<typename T>
    IndexList Indices(const BasicLookupData::Field& field,
                      const T& value,
                      const Compare::Type& compareType = Compare::EQUAL) const;

    /// Helper method that dispatches a filter request (multiple value) to the proper data member,
    /// respecting compareType.
    ///
    /// Client code (e.g. PBI filters) should use this when possible, only touching
    /// raw fields for more complex operations.
    ///
    template<typename T>
    IndexList IndicesMulti(const BasicLookupData::Field& field,
                           const std::vector<T>& values) const;

    /// \returns the virtual offsets for all records
    const std::vector<int64_t>& VirtualFileOffsets(void) const;

public:

    // map ordering doesn't make sense, optimize for direct lookup
    UnorderedLookup<int32_t> rgId_;

    // numeric comparisons make sense, keep key ordering preserved
    OrderedLookup<int32_t>  qStart_;
    OrderedLookup<int32_t>  qEnd_;
    OrderedLookup<int32_t>  holeNumber_;
    OrderedLookup<float>    readQual_;

    // see if this works, or if can use unordered, 'direct' query
    OrderedLookup<uint8_t> ctxtFlag_;

    // offsets
    std::vector<int64_t> fileOffset_;
};

// -----------------
// Mapped Data
// -----------------

class PBBAM_EXPORT MappedLookupData
{
public:
    enum Field
    {
        T_ID
      , T_START
      , T_END
      , A_START
      , A_END
      , N_M
      , N_MM
      , N_INS
      , N_DEL
      , MAP_QUALITY
      , STRAND
    };

public:
    MappedLookupData(void);
    MappedLookupData(const PbiRawMappedData& rawData);

public:

    /// Helper method that dispatches a filter request (multiple value) to the proper data member,
    /// respecting compareType.
    ///
    /// Client code (e.g. PBI filters) should use this when possible, only touching
    /// raw fields for more complex operations.
    ///
    template<typename T>
    IndexList Indices(const MappedLookupData::Field& field,
                      const T& value,
                      const Compare::Type& compareType = Compare::EQUAL) const;

    /// Helper method that dispatches a filter request (multiple value) to the proper data member,
    /// respecting compareType.
    ///
    /// Client code (e.g. PBI filters) should use this when possible, only touching
    /// raw fields for more complex operations.
    ///
    template<typename T>
    IndexList IndicesMulti(const MappedLookupData::Field& field,
                           const std::vector<T>& values) const;

public:

    // numeric comparisons make sense, keep key ordering preserved
    OrderedLookup<int32_t>  tId_;
    OrderedLookup<uint32_t> tStart_;
    OrderedLookup<uint32_t> tEnd_;
    OrderedLookup<uint32_t> aStart_;
    OrderedLookup<uint32_t> aEnd_;
    OrderedLookup<uint32_t> nM_;
    OrderedLookup<uint32_t> nMM_;
    OrderedLookup<uint8_t>  mapQV_;

    // generated values, not stored directly in PBI file
    OrderedLookup<uint32_t> nIns_;
    OrderedLookup<uint32_t> nDel_;

    // no need for map overhead, just store direct indices
    IndexList reverseStrand_;
    IndexList forwardStrand_;
};

// ------------------
// Reference Data
// ------------------

class PBBAM_EXPORT ReferenceLookupData
{
public:
    ReferenceLookupData(void);
    ReferenceLookupData(const PbiRawReferenceData& rawData);

public:

    /// \returns the range of indices (begin, end) for records that map to \p tId
    ///
    /// Client code (e.g. PBI filters) should use this when possible, only touching
    /// raw data for more complex operations.
    ///
    IndexRange Indices(const int32_t tId) const;

public:
    // references_[tId] = (begin, end) indices
    // into SubreadLookupData::fileOffset_
    std::unordered_map<int32_t, IndexRange> references_;
};

// ---------------
// Barcode Data
// ---------------

class PBBAM_EXPORT BarcodeLookupData
{
public:
    enum Field
    {
        BC_FORWARD
      , BC_REVERSE
      , BC_QUALITY
    };

public:
    // ctors
    BarcodeLookupData(void);
    BarcodeLookupData(const PbiRawBarcodeData& rawData);

public:

    /// Helper method that dispatches a filter request (multiple value) to the proper data member,
    /// respecting compareType.
    ///
    /// Client code (e.g. PBI filters) should use this when possible, only touching
    /// raw fields for more complex operations.
    ///
    template<typename T>
    IndexList Indices(const BarcodeLookupData::Field& field,
                      const T& value,
                      const Compare::Type& compareType = Compare::EQUAL) const;

    /// Helper method that dispatches a filter request (multiple value) to the proper data member,
    /// respecting compareType.
    ///
    /// Client code (e.g. PBI filters) should use this when possible, only touching
    /// raw fields for more complex operations.
    ///
    template<typename T>
    IndexList IndicesMulti(const BarcodeLookupData::Field& field,
                           const std::vector<T>& values) const;

public:
    // numeric comparisons make sense, keep key ordering preserved
    OrderedLookup<int16_t> bcForward_;
    OrderedLookup<int16_t> bcReverse_;
    OrderedLookup<int8_t>  bcQual_;
};

} // namespace BAM
} // namespace PacBio

#include "internal/PbiLookupData.inl"

#endif // PBILOOKUPDATA_H
