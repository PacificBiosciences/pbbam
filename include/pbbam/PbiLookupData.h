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
/// \file PbiLookupData.h
/// \brief Defines the classes used for PBI data lookup.
//
// Author: Derek Barnett

#ifndef PBILOOKUPDATA_H
#define PBILOOKUPDATA_H

#include "pbbam/Config.h"
#include "pbbam/Compare.h"
#include "pbbam/PbiBasicTypes.h"
#include <cstddef>
#include <cstdint>
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

/// \brief The OrderedLookup class provides a quick lookup structure for
///        PBI index data, where key values are sorted.
///
/// The main, underlying lookup structure is essentailly a std::map, where the
/// key is some value (e.g. readAccuracy) and the value is the list of indices
/// (i-th record) in the %BAM file.
///
/// This lookup class is one of the main building blocks for the PBI index
/// lookup components.
///
/// \param T    type of key stored (Accuracy for readAccuracy, int32_t for ZMW,
///             etc.)
///
template<typename T>
class OrderedLookup
{
public:
    using key_type       = T;
    using value_type     = IndexList;
    using container_type = std::map<key_type, value_type>;
    using iterator       = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty OrderedLookup structure.
    ///
    OrderedLookup() = default;

    /// \brief Creates an OrderedLookup struture, from another's underlying
    ///        lookup container.
    ///
    /// \param[in] data     lookup data container
    ///
    OrderedLookup(const container_type& data);

    /// \brief Creates an OrderedLookup struture, from another's underlying
    ///        lookup container.
    ///
    /// \param[in] data     lookup data container
    ///
    OrderedLookup(container_type&& data);

    /// \brief Creates an OrderedLookup struture, from raw data.
    ///
    /// \param[in] rawData  raw data values, where i is the index into the %BAM
    ///                     file, and rawData[i] is the key value
    ///
    OrderedLookup(const std::vector<T>& rawData);

    /// \brief Creates an OrderedLookup struture, from raw data.
    ///
    /// \param[in] rawData  raw data values, where i is the index into the %BAM
    ///                     file, and rawData[i] is the key value
    ///
    OrderedLookup(std::vector<T>&& rawData);

    /// \}

public:
    /// \name Operators
    /// \{

    /// \returns true if this lookup is same as \p other
    bool operator==(const OrderedLookup<T>& other) const;

    /// \returns true if this lookup is not the same as \p other
    bool operator!=(const OrderedLookup<T>& other) const;

    /// \}

public:
    /// \name STL-Compatibility Methods
    /// \{

    /// \returns an iterator to the first element in the underlying container
    iterator begin();

    /// \returns a const iterator to the first element in the underlying
    ///          container
    const_iterator begin() const;

    /// \returns a const iterator to the first element in the underlying
    ///
    const_iterator cbegin() const;

    /// \returns an iterator after the last element in the underlying container
    iterator end();

    /// \returns a const iterator after the last element in the underlying
    ///          container
    const_iterator end() const;

    /// \returns a const iterator after the last element in the underlying
    ///          container
    const_iterator cend() const;

    /// \returns true if underlying container is empty
    bool empty() const;

    /// \returns number of keys in the container
    size_t size() const;

    /// \}

public:
    /// \name Lookup Data
    /// \{

    /// \brief Performs a lookup into the underlying data.
    ///
    /// \param[in] key      key value to lookup
    /// \param[in] compare  compare type
    ///
    /// \returns sorted list of unique indices that satisfy the lookup key &
    ///          compare type
    ///
    IndexList LookupIndices(const key_type& key,
                            const Compare::Type& compare) const;

    /// \brief Converts the lookup structure back into its raw data.
    ///
    /// \returns raw data values, where i is the index into the %BAM file, and
    ///          rawData[i] is the key value
    ///
    std::vector<T> Unpack() const;

    /// \}

private:
    IndexList LookupInclusiveRange(const const_iterator& begin,
                                   const const_iterator& end) const;

    IndexList LookupExclusiveRange(const const_iterator& begin,
                                   const const_iterator& end,
                                   const key_type& key) const;

private:
    container_type data_;
};

/// \brief The UnorderedLookup class provides a quick lookup structure for
///        PBI index data, where key values are not sorted.
///
/// The main, underlying lookup structure is essentailly a std::unordered_map,
/// where the key is some value (e.g. read group ID) and the value is the list
/// of indices (i-th record) in the %BAM file.
///
/// This lookup class is one of the main building blocks for the PBI index
/// lookup components.
///
/// \param T    type of key stored (Accuracy for readAccuracy, int32_t for ZMW,
///             etc.)
///
template<typename T>
class UnorderedLookup
{
public:
    using key_type       = T;
    using value_type     = IndexList;
    using container_type = std::unordered_map<key_type, value_type>;
    using iterator       = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty UnorderedLookup structure.
    ///
    UnorderedLookup() = default;

    /// \brief Creates an UnorderedLookup struture, from another's underlying
    ///        lookup container.
    ///
    /// \param[in] data     lookup data container
    ///
    UnorderedLookup(const container_type& data);

    /// \brief Creates an UnorderedLookup struture, from another's underlying
    ///        lookup container.
    ///
    /// \param[in] data     lookup data container
    ///
    UnorderedLookup(container_type&& data);

    /// \brief Creates an UnorderedLookup struture, from raw data.
    ///
    /// \param[in] rawData  raw data values, where i is the index into the %BAM
    ///                     file, and rawData[i] is the key value
    ///
    UnorderedLookup(const std::vector<T>& rawData);

    /// \brief Creates an UnorderedLookup struture, from raw data.
    ///
    /// \param[in] rawData  raw data values, where i is the index into the %BAM
    ///                     file, and rawData[i] is the key value
    ///
    UnorderedLookup(std::vector<T>&& rawData);

    /// \}

public:
    /// \name Operators
    /// \{

    /// \returns true if this lookup is same as \p other
    bool operator==(const UnorderedLookup<T>& other) const;

    /// \returns true if this lookup is not the same as \p other
    bool operator!=(const UnorderedLookup<T>& other) const;

    /// \}

public:
    /// \name STL-Compatibility Methods
    /// \{

    /// \returns an iterator to the first element in the underlying container
    iterator begin();

    /// \returns a const iterator to the first element in the underlying
    ///          container
    const_iterator begin() const;

    /// \returns a const iterator to the first element in the underlying
    ///
    const_iterator cbegin() const;

    /// \returns an iterator after the last element in the underlying container
    iterator end();

    /// \returns a const iterator after the last element in the underlying
    ///          container
    const_iterator end() const;

    /// \returns a const iterator after the last element in the underlying
    ///          container
    const_iterator cend() const;

    /// \returns true if underlying container is empty
    bool empty() const;

    /// \returns number of keys in the container
    size_t size() const;

    /// \}

public:
    /// \name Lookup Data
    /// \{

    /// \brief Performs a lookup into the underlying data.
    ///
    /// \param[in] key      key value to lookup
    /// \param[in] compare  compare type
    ///
    /// \returns sorted list of unique indices that satisfy the lookup key &
    ///          compare type
    ///
    IndexList LookupIndices(const key_type& key,
                            const Compare::Type& compare) const;

    /// \brief Converts the lookup structure back into its raw data.
    ///
    /// \returns raw data values, where i is the index into the %BAM file, and
    ///          rawData[i] is the key value
    ///
    std::vector<T> Unpack() const;

    /// \}

private:
    template<typename Compare>
    IndexList LookupHelper(const key_type& key,
                           const Compare& cmp) const;

private:
    container_type data_;
};

/// \brief The BasicLookupData class provides quick lookup access to the
///        "BasicData" section of the PBI index.
///
class PBBAM_EXPORT BasicLookupData
{
public:
    /// \brief This enum describes the component fields of the BasicData
    ///        section.
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
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty lookup data object.
    BasicLookupData() = default;

    /// \brief Creates a lookup data object from the corresponding raw data.
    ///
    /// \param[in] rawData  raw data loaded from a PBI file
    ///
    BasicLookupData(const PbiRawBasicData& rawData);

    /// \}

public:
    /// \name Lookup Data Methods
    /// \{

    /// \brief Adds \b virtual file offset data to the index lookup result
    ///        blocks.
    ///
    /// A PBI lookup will result in a number of index lists, depending on the
    /// complexity of the PbiFilter involved. These index lists are then merged
    /// down into blocks of contiguous values, where each block describes a
    /// particular record index and the number of subsequent, contiguous reads
    /// that immediately follow it. In this manner, we need only perform seeks
    /// to the first record of each block.
    ///
    /// This method takes such blocks and annotates them with the corresponding
    /// \b virtual file offset. Subsequent %BAM readers can use this information
    /// to control file seeks.
    ///
    /// \param[in,out] blocks
    ///
    /// \throws std::out_of_range if a block has an invalid index value
    ///
    void ApplyOffsets(IndexResultBlocks& blocks) const;

    /// \brief This method dispatches a single-value lookup query to the proper
    ///         data member.
    ///
    /// Client code, such as custom filters, should use this when possible, only
    /// touching the raw fields for more complex operations (e.g. when unpacking
    /// is necessary).
    ///
    /// \param[in] field            section field to lookup
    /// \param[in] value            value to lookup
    /// \param[in] compareType      compare type
    ///
    /// \returns sorted list of unique indices that satisfy the lookup
    ///
    template<typename T>
    IndexList Indices(const BasicLookupData::Field& field,
                      const T& value,
                      const Compare::Type& compareType = Compare::EQUAL) const;

    /// \brief This method dispatches a multi-value lookup query to the proper
    ///        data member.
    ///
    /// Client code, such as custom filters, should use this when possible, only
    /// touching the raw fields for more complex operations (e.g. when unpacking
    /// is necessary).
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Results will correspond to an exact match on at
    ///       least one value in the list.
    ///
    /// \param[in] field        section field to lookup
    /// \param[in] values       values to lookup
    ///
    /// \returns sorted list of unique indices that satisfy the lookup
    ///
    template<typename T>
    IndexList IndicesMulti(const BasicLookupData::Field& field,
                           const std::vector<T>& values) const;

    /// \returns the \b virtual file offsets for all records
    ///
    const std::vector<int64_t>& VirtualFileOffsets() const;

    /// \}

public:
    /// \brief Lookup Data Members
    /// \{

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

    /// \}
};

/// \brief The MappedLookupData class provides quick lookup access to the
///        "MappedData" section of the PBI index.
///
class PBBAM_EXPORT MappedLookupData
{
public:
    /// \brief This enum describes the component fields of the MappedData
    ///        section.
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
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty lookup data object.
    MappedLookupData() = default;

    /// \brief Creates a lookup data object from the corresponding raw data.
    ///
    /// \param[in] rawData  raw data loaded from a PBI file
    ///
    MappedLookupData(const PbiRawMappedData& rawData);

    /// \}

public:
    /// \name Lookup Data Methods
    /// \{

    /// \brief This method dispatches a single-value lookup query to the proper
    ///         data member.
    ///
    /// Client code, such as custom filters, should use this when possible, only
    /// touching the raw fields for more complex operations (e.g. when unpacking
    /// is necessary).
    ///
    /// \param[in] field            section field to lookup
    /// \param[in] value            value to lookup
    /// \param[in] compareType      compare type
    ///
    /// \returns sorted list of unique indices that satisfy the lookup
    ///
    template<typename T>
    IndexList Indices(const MappedLookupData::Field& field,
                      const T& value,
                      const Compare::Type& compareType = Compare::EQUAL) const;

    /// \brief This method dispatches a multi-value lookup query to the proper
    ///        data member.
    ///
    /// Client code, such as custom filters, should use this when possible, only
    /// touching the raw fields for more complex operations (e.g. when unpacking
    /// is necessary).
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Results will correspond to an exact match on at
    ///       least one value in the list.
    ///
    /// \param[in] field        section field to lookup
    /// \param[in] values       values to lookup
    ///
    /// \returns sorted list of unique indices that satisfy the lookup
    ///
    template<typename T>
    IndexList IndicesMulti(const MappedLookupData::Field& field,
                           const std::vector<T>& values) const;

    /// \}

public:
    /// \name Lookup Data Members
    /// \{

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

    /// \}
};

/// \brief The ReferenceLookupData class provides quick lookup access to the
///        "CoordinateSortedData" section of the PBI index.
///
class PBBAM_EXPORT ReferenceLookupData
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty lookup data object.
    ///
    ReferenceLookupData() = default;

    /// \brief Creates a lookup data object from the corresponding raw data.
    ///
    /// \param[in] rawData  raw data loaded from a PBI file
    ///
    ReferenceLookupData(const PbiRawReferenceData& rawData);

    /// \}

public:
    /// \name Lookup Data Methods
    /// \{

    /// \brief Retrieves the index range for all records that map to a
    ///        particular reference.
    ///
    /// Client code, such as custom filters, should use this when possible, only
    /// touching the raw fields for more complex operations (e.g. when unpacking
    /// is necessary).
    ///
    /// \param[in] tId      reference ID to lookup
    ///
    /// \returns resulting index range [begin, end). If \p tId is unknown,
    ///          will return IndexRange(-1,-1) .
    ///
    IndexRange Indices(const int32_t tId) const;

    /// \}

public:
    /// \name Lookup Data Members
    /// \{

    // references_[tId] = [begin, end) indices
    std::unordered_map<int32_t, IndexRange> references_;

    /// \}
};

/// \brief The BarcodeLookupData class provides quick lookup access to the
///        "BarcodeData" section of the PBI index.
///
class PBBAM_EXPORT BarcodeLookupData
{
public:
    /// \brief This enum describes the component fields of the BarcodeData
    ///        section.
    enum Field
    {
        BC_FORWARD
      , BC_REVERSE
      , BC_QUALITY
    };

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty lookup data object.
    ///
    BarcodeLookupData() = default;

    /// \brief Creates a lookup data object from the corresponding raw data.
    ///
    /// \param[in] rawData  raw data loaded from a PBI file
    ///
    BarcodeLookupData(const PbiRawBarcodeData& rawData);

    /// \}

public:
    /// \name Lookup Data Methods
    /// \{

    /// \brief This method dispatches a single-value lookup query to the proper
    ///         data member.
    ///
    /// Client code, such as custom filters, should use this when possible, only
    /// touching the raw fields for more complex operations (e.g. when unpacking
    /// is necessary).
    ///
    /// \param[in] field            section field to lookup
    /// \param[in] value            value to lookup
    /// \param[in] compareType      compare type
    ///
    /// \returns sorted list of unique indices that satisfy the lookup
    ///
    template<typename T>
    IndexList Indices(const BarcodeLookupData::Field& field,
                      const T& value,
                      const Compare::Type& compareType = Compare::EQUAL) const;

    /// \brief This method dispatches a multi-value lookup query to the proper
    ///        data member.
    ///
    /// Client code, such as custom filters, should use this when possible, only
    /// touching the raw fields for more complex operations (e.g. when unpacking
    /// is necessary).
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Results will correspond to an exact match on at
    ///       least one value in the list.
    ///
    /// \param[in] field        section field to lookup
    /// \param[in] values       values to lookup
    ///
    /// \returns sorted list of unique indices that satisfy the lookup
    ///
    template<typename T>
    IndexList IndicesMulti(const BarcodeLookupData::Field& field,
                           const std::vector<T>& values) const;

    /// \}

public:
    /// \name Lookup Data Members
    /// \{

    // numeric comparisons make sense, keep key ordering preserved
    OrderedLookup<int16_t> bcForward_;
    OrderedLookup<int16_t> bcReverse_;
    OrderedLookup<int8_t>  bcQual_;

    /// \}
};

} // namespace BAM
} // namespace PacBio

#include "internal/PbiLookupData.inl"

#endif // PBILOOKUPDATA_H
