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
/// \file PbiFilterTypes.h
/// \brief Defines the built-in PBI filters.
//
// Author: Derek Barnett

#ifndef PBIFILTERTYPES_H
#define PBIFILTERTYPES_H

#include "pbbam/Compare.h"
#include "pbbam/PbiFilter.h"
#include "pbbam/PbiIndex.h"
#include <boost/optional.hpp>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

namespace internal {

/// \internal
///
/// Provides basic container for value/compare-type pair
///
template<typename T>
struct FilterBase
{
public:
    T value_;
    boost::optional<std::vector<T> > multiValue_;
    Compare::Type cmp_;
protected:
    FilterBase(const T& value, const Compare::Type cmp);
    FilterBase(T&& value, const Compare::Type cmp);
    FilterBase(const std::vector<T>& values);
    FilterBase(std::vector<T>&& values);
protected:
    bool CompareHelper(const T& lhs) const;
private:
    bool CompareSingleHelper(const T& lhs) const;
    bool CompareMultiHelper(const T& lhs) const;
};

/// \internal
///
/// Dispatches the lookup to BarcodeLookupData
///
template<typename T, BarcodeLookupData::Field field>
struct BarcodeDataFilterBase : public FilterBase<T>
{
protected:
    BarcodeDataFilterBase(const T& value, const Compare::Type cmp);
    BarcodeDataFilterBase(T&& value, const Compare::Type cmp);
    BarcodeDataFilterBase(const std::vector<T>& values);
    BarcodeDataFilterBase(std::vector<T>&& values);
public:
    bool Accepts(const PbiRawData& idx, const size_t row) const;
};

/// \internal
///
/// Dispatches the lookup to BasicLookupData
///
template<typename T, BasicLookupData::Field field>
struct BasicDataFilterBase : public FilterBase<T>
{
protected:
    BasicDataFilterBase(const T& value, const Compare::Type cmp);
    BasicDataFilterBase(T&& value, const Compare::Type cmp);
    BasicDataFilterBase(const std::vector<T>& values);
    BasicDataFilterBase(std::vector<T>&& values);
public:
    bool Accepts(const PbiRawData& idx, const size_t row) const;
};

/// \internal
///
/// Dispatches the lookup to MappedLookupData
///
template<typename T, MappedLookupData::Field field>
struct MappedDataFilterBase : public FilterBase<T>
{
protected:
    MappedDataFilterBase(const T& value, const Compare::Type cmp);
    MappedDataFilterBase(T&& value, const Compare::Type cmp);
    MappedDataFilterBase(const std::vector<T>& values);
    MappedDataFilterBase(std::vector<T>&& values);
public:
    bool Accepts(const PbiRawData& idx, const size_t row) const;
};

} // namespace internal

/// \brief The PbiAlignedEndFilter class provides a PbiFilter-compatible filter
///        on aligned end.
///
/// Example: \include code/PbiAlignedEndFilter.txt
///
/// \sa BamRecord::AlignedEnd
///
struct PbiAlignedEndFilter
    : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::A_END>
{
public:
    /// \brief Creates a filter on aligned end.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiAlignedEndFilter(const uint32_t position,
                        const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiAlignedLengthFilter class provides a PbiFilter-compatible
///        filter on aligned length.
///
/// Example: \include code/PbiAlignedLengthFilter.txt
///
/// \sa BamRecord::AlignedEnd, BamRecord::AlignedStart
///
struct PbiAlignedLengthFilter : public internal::FilterBase<uint32_t>
{
public:
    /// \brief Creates a filter on aligned length.
    ///
    /// \param[in] length value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiAlignedLengthFilter(const uint32_t length,
                           const Compare::Type cmp = Compare::EQUAL);

public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;
};

/// \brief The PbiAlignedStartFilter class provides a PbiFilter-compatible
///        filter on aligned start.
///
/// Example: \include code/PbiAlignedStartFilter.txt
///
/// \sa BamRecord::AlignedStart
///
struct PbiAlignedStartFilter
    : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::A_START>
{
public:
    /// \brief Creates a filter on aligned start.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiAlignedStartFilter(const uint32_t position,
                          const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiAlignedStrandFilter class provides a PbiFilter-compatible
///        filter on aligned strand.
///
/// Example: \include code/PbiAlignedStrandFilter.txt
///
/// \sa BamRecord::AlignedStrand
///
struct PbiAlignedStrandFilter
    : public internal::MappedDataFilterBase<Strand, MappedLookupData::STRAND>
{
public:
    /// \brief Creates a strand filter.
    ///
    /// \param[in] strand  strand value to compare on
    /// \param[in] cmp     compare type
    ///
    PbiAlignedStrandFilter(const Strand strand,
                           const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiBarcodeFilter class provides a PbiFilter-compatible filter on
///        barcode ID.
///
/// Any record with this barcode ID (forward or reverse) will pass this filter.
///
/// Example: \include code/PbiBarcodeFilter.txt
///
/// \sa BamRecord::BarcodeForward, BamRecord::BarcodeReverse
///
struct PbiBarcodeFilter
{
public:
    /// \brief Creates a single-value barcode filter.
    ///
    /// \param[in] barcode  barcode ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeFilter(const int16_t barcode,
                     const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a 'whitelisted' barcode filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly, in either bc_forward or bc_reverse.
    ///
    /// \param[in] whitelist  barcode IDs to compare on
    ///
    PbiBarcodeFilter(const std::vector<int16_t>& whitelist);

    /// \brief Creates a 'whitelisted' barcode filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly, in either bc_forward or bc_reverse.
    ///
    /// \param[in] whitelist  barcode IDs to compare on
    ///
    PbiBarcodeFilter(std::vector<int16_t>&& whitelist);

public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    PbiFilter compositeFilter_;
};

/// \brief The PbiBarcodeForwardFilter class provides a PbiFilter-compatible
///        filter on forward barcode ID.
///
/// Example: \include code/PbiBarcodeForwardFilter.txt
///
/// \sa BamRecord::BarcodeForward
///
struct PbiBarcodeForwardFilter
    : public internal::BarcodeDataFilterBase<int16_t, BarcodeLookupData::BC_FORWARD>
{
public:
    /// \brief Creates a single-value forward barcode filter.
    ///
    /// \param[in] bcFwdId  (forward) barcode ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeForwardFilter(const int16_t bcFwdId,
                            const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a 'whitelisted' forward barcode filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly, in bc_forward.
    ///
    /// \param[in] whitelist  barcode IDs to compare on
    ///
    PbiBarcodeForwardFilter(const std::vector<int16_t>& whitelist);

    /// \brief Creates a 'whitelisted' forward barcode filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly, in bc_forward.
    ///
    /// \param[in] whitelist  barcode IDs to compare on
    ///
    PbiBarcodeForwardFilter(std::vector<int16_t>&& whitelist);
};

/// \brief The PbiBarcodeQualityFilter class provides a PbiFilter-compatible
///        filter on  barcode quality.
///
/// Example: \include code/PbiBarcodeQualityFilter.txt
///
/// \sa BamRecord::BarcodeQuality
///
struct PbiBarcodeQualityFilter
    : public internal::BarcodeDataFilterBase<uint8_t, BarcodeLookupData::BC_QUALITY>
{
public:
    /// \brief Creates a single-value barcode quality filter.
    ///
    /// \param[in] bcQuality    barcode quality to compare on
    /// \param[in] cmp          compare type
    ///
    PbiBarcodeQualityFilter(const uint8_t bcQuality,
                            const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiBarcodeReverseFilter class provides a PbiFilter-compatible
///        filter on forward barcode ID.
///
/// Example: \include code/PbiBarcodeReverseFilter.txt
///
/// \sa BamRecord::BarcodeReverse
///
struct PbiBarcodeReverseFilter
    : public internal::BarcodeDataFilterBase<int16_t, BarcodeLookupData::BC_REVERSE>
{
public:
    /// \brief Creates a single-value reverse barcode filter.
    ///
    /// \param[in] bcRevId  (reverse) barcode ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeReverseFilter(const int16_t bcRevId,
                            const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a 'whitelisted' reverse barcode filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly, in bc_reverse.
    ///
    /// \param[in] whitelist  barcode IDs to compare on
    ///
    PbiBarcodeReverseFilter(const std::vector<int16_t>& whitelist);

    /// \brief Creates a 'whitelisted' reverse barcode filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly, in bc_reverse.
    ///
    /// \param[in] whitelist  barcode IDs to compare on
    ///
    PbiBarcodeReverseFilter(std::vector<int16_t>&& whitelist);
};

/// \brief The PbiBarcodesFilter class provides a PbiFilter-compatible filter on
///        both forward & reverse barcode IDs.
///
/// A record must match both IDs to pass the filter.
///
/// Example: \include code/PbiBarcodesFilter.txt
///
/// \sa BamRecord::Barcodes
///
struct PbiBarcodesFilter
{
public:
    /// \brief Creates a barcodes filter from a std::pair of IDs.
    ///
    /// pair.first -> BarcodeForward\n
    /// pair.second -> BarcodeReverse
    ///
    /// \param[in] barcodes barcode IDs to compare on
    /// \param[in] cmp      compare type
    ///
    PbiBarcodesFilter(const std::pair<int16_t, int16_t> barcodes,
                      const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a barcodes filter from forward & reverse IDs.
    ///
    /// \param[in] bcForward    forward barcode ID to compare on
    /// \param[in] bcReverse    reverse barcode ID to compare on
    /// \param[in] cmp          compare type
    ///
    PbiBarcodesFilter(const int16_t bcForward,
                      const int16_t bcReverse,
                      const Compare::Type cmp = Compare::EQUAL);
public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    PbiFilter compositeFilter_;
};

/// \brief The PbiIdentityFilter class provides a PbiFilter-compatible filter on
///        read identity (% aligned match).
///
/// Read identity is equivalent to: 1.0 - (nMM + nDel + nIns)/readLength.
///
/// Example: \include code/PbiIdentityFilter.txt
///
struct PbiIdentityFilter : public internal::FilterBase<float>
{
public:
    /// \brief Creates a read identity filter.
    ///
    /// \param[in] identity value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiIdentityFilter(const float identity,
                      const Compare::Type cmp = Compare::EQUAL);

public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;
};

/// \brief The PbiLocalContextFilter class provides a PbiFilter-compatible
///        filter on local context (adapter, barcode, etc.).
///
/// The primary Compare::Type operators intended for this filter are:
/// Compare::EQUAL, Compare::NOT_EQUAL, Compare::CONTAINS, and
/// Compare::NOT_CONTAINS.
///
/// Example: \include code/PbiLocalContextFilter.txt
///
struct PbiLocalContextFilter
    : public internal::BasicDataFilterBase<LocalContextFlags,
                                           BasicLookupData::CONTEXT_FLAG >
{
public:
    PbiLocalContextFilter(const LocalContextFlags& flags,
                          const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiMapQualityFilter class provides a PbiFilter-compatible filter on
///        mapping quality.
///
/// Example: \include code/PbiMapQualityFilter.txt
///
/// \sa BamRecord::MapQuality
///
struct PbiMapQualityFilter
    : public internal::MappedDataFilterBase<uint8_t, MappedLookupData::MAP_QUALITY>
{
public:
    /// \brief Creates a map quality filter.
    ///
    /// \param[in] mapQual  value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiMapQualityFilter(const uint8_t mapQual,
                        const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiMovieNameFilter class provides a PbiFilter-compatible filter
///        on movie name.
///
/// Example: \include code/PbiMovieNameFilter.txt
///
/// \sa BamRecord::MovieName
///
struct PbiMovieNameFilter
{
public:
    /// \brief Creates a single-value movie name filter.
    ///
    /// \param[in] movieName    movie name to compare on
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match movie name, exactly.
    ///
    PbiMovieNameFilter(const std::string& movieName);

    /// \brief Creates a 'whitelisted' movie name filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    movie names to compare on
    ///
    PbiMovieNameFilter(const std::vector<std::string>& whitelist);

    /// \brief Creates a 'whitelisted' movie name filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    movie names to compare on
    ///
    PbiMovieNameFilter(std::vector<std::string>&& whitelist);

public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
   PbiFilter compositeFilter_;
};

/// \brief The PbiNumDeletedBasesFilter class provides a PbiFilter-compatible
///        filter on the number of deleted bases.
///
/// Example: \include code/PbiNumDeletedBasesFilter.txt
///
/// \sa BamRecord::NumDeletedBases
///
struct PbiNumDeletedBasesFilter
    : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_DEL>
{
public:
    /// \brief Creates a filter on the number of deleted bases.
    ///
    /// \param[in] numDeletions value to compare on
    /// \param[in] cmp          compare type
    ///
    PbiNumDeletedBasesFilter(const size_t numDeletions,
                             const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiNumInsertededBasesFilter class provides a PbiFilter-compatible
///        filter on the number of inserted bases.
///
/// Example: \include code/PbiNumInsertedBasesFilter.txt
///
/// \sa BamRecord::NumInsertedBases
///
struct PbiNumInsertedBasesFilter
    : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_INS>
{
public:
    /// \brief Creates a filter on the number of inserted bases.
    ///
    /// \param[in] numInsertions    value to compare on
    /// \param[in] cmp              compare type
    ///
    PbiNumInsertedBasesFilter(const size_t numInsertions,
                              const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiNumMatchesFilter class provides a PbiFilter-compatible filter
///        on the number of matched bases.
///
/// Example: \include code/PbiNumMatchesFilter.txt
///
/// \sa BamRecord::NumMatches
///
struct PbiNumMatchesFilter
    : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_M>
{
public:
    /// \brief Creates a filter on the number of matched bases.
    ///
    /// \param[in] numMatchedBases  value to compare on
    /// \param[in] cmp              compare type
    ///
    PbiNumMatchesFilter(const size_t numMatchedBases,
                        const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiNumMismatchesFilter class provides a PbiFilter-compatible
///        filter on the number of mismatched bases.
///
/// Example: \include code/PbiNumMismatchesFilter.txt
///
/// \sa BamRecord::NumMismatches
///
struct PbiNumMismatchesFilter
    : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_MM>
{
public:
    /// \brief Creates a filter on the number of mismatched bases.
    ///
    /// \param[in] numMismatchedBases   value to compare on
    /// \param[in] cmp                  compare type
    ///
    PbiNumMismatchesFilter(const size_t numMismatchedBases,
                           const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiQueryEndFilter class provides a PbiFilter-compatible filter
///        on query end.
///
/// Example: \include code/PbiQueryEndFilter.txt
///
/// \sa BamRecord::QueryEnd
///
struct PbiQueryEndFilter
    : public internal::BasicDataFilterBase<int32_t, BasicLookupData::Q_END>
{
public:
    /// \brief Creates a filter on query end position.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiQueryEndFilter(const int32_t position,
                      const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiQueryLengthFilter class provides a PbiFilter-compatible filter
///        on query length.
///
/// queryLength = (queryEnd - queryStart)
///
/// Example: \include code/PbiQueryLengthFilter.txt
///
/// \sa BamRecord::QueryEnd, BamRecord::QueryStart
///
struct PbiQueryLengthFilter : public internal::FilterBase<int32_t>
{
public:
    /// \brief Creates a filter on query length
    ///
    /// \param[in] length   value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiQueryLengthFilter(const int32_t length,
                         const Compare::Type cmp = Compare::EQUAL);

public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;
};

/// \brief The PbiQueryNameFilter class provides a PbiFilter-compatible filter
///        on name length.
///
/// Example: \include code/PbiQueryNameFilter.txt
///
/// \sa BamRecord::FullName
///
struct PbiQueryNameFilter
{
public:
    /// \brief Creates a single-value query name filter.
    ///
    /// \param[in] qname    query name to compare on
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match query name, exactly.
    ///
    PbiQueryNameFilter(const std::string& qname);

    /// \brief Creates a 'whitelisted' query name filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    query names to compare on
    ///
    PbiQueryNameFilter(const std::vector<std::string>& whitelist);

    PbiQueryNameFilter(const PbiQueryNameFilter& other);
    ~PbiQueryNameFilter();

public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    struct PbiQueryNameFilterPrivate;
    std::unique_ptr<PbiQueryNameFilterPrivate> d_;
};

/// \brief The PbiQueryStartFilter class provides a PbiFilter-compatible filter
///        on query start.
///
/// Example: \include code/PbiQueryStartFilter.txt
///
/// \sa BamRecord::QueryStart
///
struct PbiQueryStartFilter
    : public internal::BasicDataFilterBase<int32_t, BasicLookupData::Q_START>
{
public:
    /// \brief Creates a filter on query start position.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiQueryStartFilter(const int32_t position,
                        const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiReadAccuracyFilter class provides a PbiFilter-compatible filter
///        on read accuracy.
///
/// Example: \include code/PbiReadAccuracyFilter.txt
///
/// \sa BamRecord::ReadAccuracy
///
struct PbiReadAccuracyFilter
    : public internal::BasicDataFilterBase<Accuracy, BasicLookupData::READ_QUALITY>
{
public:
    /// \brief Creates a filter on read accuracy.
    ///
    /// \param[in] accuracy value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReadAccuracyFilter(const Accuracy accuracy,
                          const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiReadGroupFilter class provides a PbiFilter-compatible filter
///        on read group.
///
/// Example: \include code/PbiReadGroupFilter.txt
///
/// \sa BamRecord::ReadGroup,
///     BamRecord::ReadGroupId,
///     BamRecord::ReadGroupNumericId
///
struct PbiReadGroupFilter
    : public internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>
{
public:
    /// \brief Creates a filter on read group (numeric) ID value
    ///
    /// \param[in] rgId     numeric read group ID
    /// \param[in] cmp      compare type
    ///
    /// \sa BamRecord::ReadGroupNumericId
    ///
    PbiReadGroupFilter(const int32_t rgId,
                       const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a filter on printable read group ID value
    ///
    /// \param[in] rgId     read group ID string
    /// \param[in] cmp      compare type
    ///
    /// \sa BamRecord::ReadGroupId
    ///
    PbiReadGroupFilter(const std::string rgId,
                       const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a filter on read group (object).
    ///
    /// \param[in] rg   read group object
    /// \param[in] cmp  compare type
    ///
    /// \sa BamRecord::ReadGroup
    ///
    PbiReadGroupFilter(const ReadGroupInfo& rg,
                       const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a 'whitelisted' filter on read group numeric IDs.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    read group IDs to compare on
    ///
    PbiReadGroupFilter(const std::vector<int32_t>& whitelist);

    /// \brief Creates a 'whitelisted' filter on read group numeric IDs.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    read group IDs to compare on
    ///
    PbiReadGroupFilter(std::vector<int32_t>&& whitelist);

    /// \brief Creates a 'whitelisted' filter on read group printable IDs.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    read group ID strings to compare on
    ///
    PbiReadGroupFilter(const std::vector<std::string>& whitelist);

    /// \brief Creates a 'whitelisted' filter on read group printable IDs.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    read group ID strings to compare on
    ///
    PbiReadGroupFilter(std::vector<std::string>&& whitelist);

    /// \brief Creates a 'whitelisted' filter using read group objects.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    read group objects to compare on
    ///
    PbiReadGroupFilter(const std::vector<ReadGroupInfo>& whitelist);

    /// \brief Creates a 'whitelisted' filter using read group objects.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    read group objects to compare on
    ///
    PbiReadGroupFilter(std::vector<ReadGroupInfo>&& whitelist);
};

/// \brief The PbiReferenceEndFilter class provides a PbiFilter-compatible
///        filter on reference end.
///
/// Example: \include code/PbiReferenceEndFilter.txt
///
/// \sa BamRecord::ReferenceEnd
///
struct PbiReferenceEndFilter
    : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::T_END>
{
public:
    /// \brief Creates a filter on reference end.
    ///
    /// \param[in] tEnd     value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReferenceEndFilter(const uint32_t tEnd,
                          const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiReferenceIdFilter class provides a PbiFilter-compatible
///        filter on reference ID.
///
/// Example: \include code/PbiReferenceIdFilter.txt
///
/// \sa BamRecord::ReferenceId
///
struct PbiReferenceIdFilter
    : public internal::MappedDataFilterBase<int32_t, MappedLookupData::T_ID>
{
public:
    /// \brief Creates a single-value reference ID filter.
    ///
    /// \param[in] tId  reference ID to compare on
    /// \param[in] cmp  compare type
    ///
    PbiReferenceIdFilter(const int32_t tId,
                         const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a 'whitelisted' reference ID filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    reference IDs to compare on
    ///
    PbiReferenceIdFilter(const std::vector<int32_t>& whitelist);

    /// \brief Creates a 'whitelisted' reference ID filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    reference IDs to compare on
    ///
    PbiReferenceIdFilter(std::vector<int32_t>&& whitelist);
};

/// \brief The PbiReferenceNameFilter class provides a PbiFilter-compatible
///        filter on reference name.
///
/// Example: \include code/PbiReferenceNameFilter.txt
///
/// \sa BamRecord::ReferenceName
///
struct PbiReferenceNameFilter
{
public:
    /// \brief Creates a single-value reference name filter.
    ///
    /// \param[in] rname    reference ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReferenceNameFilter(std::string rname,
                           Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a 'whitelisted' reference name filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    reference names to compare on
    ///
    PbiReferenceNameFilter(const std::vector<std::string>& whitelist);

    /// \brief Creates a 'whitelisted' reference name filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    reference names to compare on
    ///
    PbiReferenceNameFilter(std::vector<std::string>&& whitelist);

public:
    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    mutable bool initialized_ = false;
    mutable PbiFilter subFilter_;
    std::string rname_;
    boost::optional<std::vector<std::string> > rnameWhitelist_;
    Compare::Type cmp_;

private:
    // marked const so we can delay setup of filter in Accepts(), once we have
    // access to PBI/BAM input. modified values marked mutable accordingly
    void Initialize(const PbiRawData& idx) const;
};

/// \brief The PbiReferenceStartFilter class provides a PbiFilter-compatible
///        filter on reference start.
///
/// Example: \include code/PbiReferenceStartFilter.txt
///
/// \sa BamRecord::ReferenceStart
///
struct PbiReferenceStartFilter
    : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::T_START>
{
public:
    /// \brief Creates a filter on reference start.
    ///
    /// \param[in] tStart   value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReferenceStartFilter(const uint32_t tStart,
                            const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiZmwFilter class provides a PbiFilter-compatible filter on
///        ZMW hole number.
///
/// Example: \include code/PbiZmwFilter.txt
///
/// \sa BamRecord::HoleNumber
///
struct PbiZmwFilter : public internal::BasicDataFilterBase<int32_t,
                                                           BasicLookupData::ZMW>
{
public:
    /// \brief Creates a single-value ZMW hole number filter.
    ///
    /// \param[in] zmw  value to compare on
    /// \param[in] cmp  compare type
    ///
    PbiZmwFilter(const int32_t zmw,
                 const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a 'whitelisted' ZMW hole number filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    ZMW hole numbers to compare on
    ///
    PbiZmwFilter(const std::vector<int32_t>& whitelist);

    /// \brief Creates a 'whitelisted' ZMW hole number filter.
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match at least one value from the
    ///       whitelist, exactly.
    ///
    /// \param[in] whitelist    ZMW hole numbers to compare on
    ///
    PbiZmwFilter(std::vector<int32_t>&& whitelist);
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/PbiFilterTypes.inl"

#endif // PBIFILTERTYPES_H
