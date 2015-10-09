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

#ifndef PBIFILTERTYPES_H
#define PBIFILTERTYPES_H

#include "pbbam/Compare.h"
#include "pbbam/PbiFilter.h"
#include "pbbam/PbiIndex.h"
#include <boost/optional.hpp>
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
protected:
    FilterBase(const T& value, const Compare::Type cmp);
    FilterBase(T&& value, const Compare::Type cmp);
    FilterBase(const std::vector<T>& values);
    FilterBase(std::vector<T>&& values);

public:
    T value_;
    boost::optional<std::vector<T> > multiValue_;
    Compare::Type cmp_;
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
    IndexList Lookup(const PbiIndex& index) const;
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
    IndexList Lookup(const PbiIndex& index) const;
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
    IndexList Lookup(const PbiIndex& index) const;
};

} // namespace internal

/// Provides a filter that operates on PBI field MappedData::aEnd
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiAlignedEndFilter{3000, Compare::GREATER_THAN} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiAlignedEndFilter : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::A_END>
{
    PbiAlignedEndFilter(const uint32_t position, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on aligned length (aEnd - aStart)
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiQueryLengthFilter{ 2000 , Compare::GREATER_THAN } };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiAlignedLengthFilter : public internal::FilterBase<uint32_t>
{
    PbiAlignedLengthFilter(const uint32_t length, const Compare::Type cmp = Compare::EQUAL);

    IndexList Lookup(const PbiIndex& idx) const;
};

/// Provides a filter that operates on PBI field MappedData::aStart
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiAlignedStartFilter{3000, Compare::GREATER_THAN} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiAlignedStartFilter : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::A_START>
{
    PbiAlignedStartFilter(const uint32_t position, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field MappedData::revStrand
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiAlignedStrandFilter{Strand::FORWARD} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiAlignedStrandFilter : public internal::MappedDataFilterBase<Strand, MappedLookupData::STRAND>
{
    PbiAlignedStrandFilter(const Strand strand, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on an barcode ID. Any record with this barcode ID (forward or reverse) will pass this filter.
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiBarcodeFilter{17} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
/// \note The whitelist versions only use Compare::EQUAL. Records will match at least one value
///       provided in the whitelist, exactly, in either bc_forward or bc_reverse.
///
struct PbiBarcodeFilter
{
public:
    PbiBarcodeFilter(const uint16_t barcode, const Compare::Type cmp = Compare::EQUAL);
    PbiBarcodeFilter(const std::vector<uint16_t>& whitelist);
    PbiBarcodeFilter(std::vector<uint16_t>&& whitelist);

    IndexList Lookup(const PbiIndex& idx) const;

private:
    PbiFilter compositeFilter_;
};

/// Provides a filter that operates on PBI field BarcodeData::bc_forward.
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiBarcodeForwardFilter{50} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
/// \note The whitelist versions only use Compare::EQUAL. Records will match at least one value
///       provided in the whitelist, exactly, in bc_forward.
///
struct PbiBarcodeForwardFilter : public internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_FORWARD>
{
    PbiBarcodeForwardFilter(const uint16_t bcFwdId, const Compare::Type cmp = Compare::EQUAL);
    PbiBarcodeForwardFilter(const std::vector<uint16_t>& whitelist);
    PbiBarcodeForwardFilter(std::vector<uint16_t>&& whitelist);
};

/// Provides a filter that operates on PBI field BarcodeData::bc_qual
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiBarcodeForwardFilter{42, Compare::GREATER_THAN_EQUAL} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiBarcodeQualityFilter : public internal::BarcodeDataFilterBase<uint8_t, BarcodeLookupData::BC_QUALITY>
{
    PbiBarcodeQualityFilter(const uint8_t bcQuality, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field BarcodeData::bc_reverse
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiBarcodeReverseFilter{50} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
/// \note The whitelist versions only use Compare::EQUAL. Records will match at least one value
///       provided in the whitelist, exactly, in bc_reverse.
///
struct PbiBarcodeReverseFilter : public internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_REVERSE>
{
    PbiBarcodeReverseFilter(const uint16_t bcRevId, const Compare::Type cmp = Compare::EQUAL);
    PbiBarcodeReverseFilter(const std::vector<uint16_t>& whitelist);
    PbiBarcodeReverseFilter(std::vector<uint16_t>&& whitelist);
};

/// Provides a filter that operates on barcode ID pairs. A record must match both IDs to pass the filter.
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiBarcodesFilter{17, 18} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiBarcodesFilter
{
public:
    /// Constructs filter from a std::pair of barcode IDs (forward, reverse).
    ///
    /// \sa BamRecord::Barcodes
    ///
    PbiBarcodesFilter(const std::pair<uint16_t, uint16_t> barcodes,
                      const Compare::Type cmp = Compare::EQUAL);

    /// Constructs filter from a pair of barcode IDs (forward, reverse).
    ///
    /// \sa BamRecord::Barcodes
    ///
    PbiBarcodesFilter(const uint16_t bcForward,
                      const uint16_t bcReverse,
                      const Compare::Type cmp = Compare::EQUAL);

    IndexList Lookup(const PbiIndex& idx) const;

private:
    PbiFilter compositeFilter_;
};

/// Provides a filter that operates on read identity (% aligned match).
///
/// Equivalent to 1.0 - (nMM + nDel + nIns)/readLength.
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiIdentityFilter{0.9, Compare::GREATER_THAN_EQUAL} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiIdentityFilter : public internal::FilterBase<float>
{
    PbiIdentityFilter(const float identity, const Compare::Type cmp = Compare::EQUAL);
    IndexList Lookup(const PbiIndex& idx) const;
};

// TODO: determine use case(s) for query - entire flag or parts?
//struct PbiLocalContextFlagFilter : public internal::BasicDataFilterBase<LocalContextFlags, BasicLookupData::CONTEXT_FLAG > { };

/// Provides a filter that operates on PBI field MappedData::mapQV
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiMapQualityFilter{75, Compare::GREATER_THAN_EQUAL} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiMapQualityFilter : public internal::MappedDataFilterBase<uint8_t, MappedLookupData::MAP_QUALITY>
{
    PbiMapQualityFilter(const uint8_t mapQual, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on movie name.
///
/// \note Unlike most other filters, Compare::EQUAL is the only supported compare type,
/// hence no optional argument for it in the constructor.
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiMovieNameFilter{ "movie_name" } };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiMovieNameFilter : public internal::FilterBase<std::string>
{
    PbiMovieNameFilter(const std::string& movieName);
    PbiMovieNameFilter(const std::vector<std::string>& whitelist);
    PbiMovieNameFilter(std::vector<std::string>&& whitelist);

    IndexList Lookup(const PbiIndex& idx) const;
};

/// Provides a filter that operates on MappedData::nDel (not stored in PBI, but calculated from other data)
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiNumDeletedBasesFilter{50, Compare::LESS_THAN} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiNumDeletedBasesFilter : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_DEL>
{
    PbiNumDeletedBasesFilter(const size_t numDeletions, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on MappedData::nIns (not stored in PBI, but calculated from other data)
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiNumDInsertedBasesFilter{50, Compare::LESS_THAN} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiNumInsertedBasesFilter : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_INS>
{
    PbiNumInsertedBasesFilter(const size_t numInsertions, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field MappedData::nM
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiNumMatchesFilter{4000, Compare::GREATER_THAN_EQUAL} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiNumMatchesFilter : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_M>
{
    PbiNumMatchesFilter(const size_t numMatchedBases, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field MappedData::nMM
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiNumMismatchesFilter{400, Compare::LESS_THAN_EQUAL} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiNumMismatchesFilter : public internal::MappedDataFilterBase<size_t, MappedLookupData::N_MM>
{
    PbiNumMismatchesFilter(const size_t numMismatchedBases, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field BasicData::qEnd
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiQueryEndFilter{3000} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiQueryEndFilter : public internal::BasicDataFilterBase<int32_t, BasicLookupData::Q_END>
{
    PbiQueryEndFilter(const int32_t position, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on query length (qEnd - qStart)
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiQueryLengthFilter{ 2000 , Compare::GREATER_THAN } };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiQueryLengthFilter : public internal::FilterBase<int32_t>
{
    PbiQueryLengthFilter(const int32_t length, const Compare::Type cmp = Compare::EQUAL);
    IndexList Lookup(const PbiIndex& idx) const;
};

/// Provides a filter that operates on query name.
///
/// \note Unlike most other filters, Compare::EQUAL is the only supported compare type,
/// hence no optional argument for it in the constructor.
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiQueryNameFilter{ "movie_name/42/100_200" } };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiQueryNameFilter : public internal::FilterBase<std::string>
{
    PbiQueryNameFilter(const std::string& qname);
    PbiQueryNameFilter(const std::vector<std::string>& whitelist);
    PbiQueryNameFilter(std::vector<std::string>&& whitelist);

    IndexList Lookup(const PbiIndex& idx) const;
};

/// Provides a filter that operates on PBI field BasicData::qStart
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiQueryEndFilter{3000} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiQueryStartFilter : public internal::BasicDataFilterBase<int32_t, BasicLookupData::Q_START>
{
    PbiQueryStartFilter(const int32_t position, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field BasicData::readQual
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiReadAccuracyFilter{0.8, Compare::GREATER_THAN_EQUAL} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiReadAccuracyFilter : public internal::BasicDataFilterBase<Accuracy, BasicLookupData::READ_QUALITY>
{
    PbiReadAccuracyFilter(const Accuracy accuracy, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field BasicData::rgId
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiReadGroupFilter{ 2458765 } };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
/// \note The whitelist versions only use Compare::EQUAL. Records will match at least one value
///       provided in the whitelist.
///
struct PbiReadGroupFilter : public internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>
{
    /// Constructs a read group filter using the raw, numeric ID value.
    PbiReadGroupFilter(const int32_t rgId, const Compare::Type cmp = Compare::EQUAL);

    /// Constructs a read group filter using the printable string used in the BAM file.
    PbiReadGroupFilter(const std::string rgId, const Compare::Type cmp = Compare::EQUAL);

    /// Constructs a read group filter using a ReadGroupInfo object from a BamHeader.
    PbiReadGroupFilter(const ReadGroupInfo& rg, const Compare::Type cmo = Compare::EQUAL);

    /// Constructs a read group filter using raw, numeric ID values.
    PbiReadGroupFilter(const std::vector<int32_t>& whitelist);

    /// Constructs a read group filter using raw, numeric ID values.
    PbiReadGroupFilter(std::vector<int32_t>&& whitelist);

    /// Constructs a read group filter using the printable strings used in the BAM file.
    PbiReadGroupFilter(const std::vector<std::string>& whitelist);

    /// Constructs a read group filter using the printable strings used in the BAM file.
    PbiReadGroupFilter(std::vector<std::string>&& whitelist);

    /// Constructs a read group filter using ReadGroupInfo objects from a BamHeader.
    PbiReadGroupFilter(const std::vector<ReadGroupInfo>& whitelist);

    /// Constructs a read group filter using ReadGroupInfo objects from a BamHeader.
    PbiReadGroupFilter(std::vector<ReadGroupInfo>&& whitelist);
};

/// Provides a filter that operates on PBI field MappedData::tEnd
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiReferenceEndFilter{2000} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiReferenceEndFilter : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::T_END>
{
    PbiReferenceEndFilter(const uint32_t tEnd, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field MappedData::tId
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiReferenceIdFilter{6} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
/// \note The whitelist versions only use Compare::EQUAL. Records will match at least one value
///       provided in the whitelist.
///
struct PbiReferenceIdFilter : public internal::MappedDataFilterBase<int32_t, MappedLookupData::T_ID>
{
    PbiReferenceIdFilter(const int32_t tId, const Compare::Type cmp = Compare::EQUAL);
    PbiReferenceIdFilter(const std::vector<int32_t>& whitelist);
    PbiReferenceIdFilter(std::vector<int32_t>&& whitelist);
};

/// Provides a filter that operates on reference name
///
/// \note Compare::EQUAL and Compare::NOT_EQUAL are the only supported compare types.
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiReferenceNameFilter{"chr1"} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
///\note The whitelist versions only use Compare::EQUAL. Records will match at least one value
///       provided in the whitelist.
///
struct PbiReferenceNameFilter : public internal::FilterBase<std::string>
{
    PbiReferenceNameFilter(const std::string& rname, const Compare::Type cmp = Compare::EQUAL);
    PbiReferenceNameFilter(const std::vector<std::string>& whitelist);
    PbiReferenceNameFilter(std::vector<std::string>&& whitelist);

    IndexList Lookup(const PbiIndex& idx) const;
};

/// Provides a filter that operates on PBI field MappedData::tStart
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiReferenceStartFilter{2000} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
struct PbiReferenceStartFilter : public internal::MappedDataFilterBase<uint32_t, MappedLookupData::T_START>
{
    PbiReferenceStartFilter(const uint32_t tStart, const Compare::Type cmp = Compare::EQUAL);
};

/// Provides a filter that operates on PBI field BasicData::holeNumber
///
/// Example:
/// \code{.cpp}
///     PbiFilter filter{ PbiZmwFilter{42} };
///     PbiFilterQuery query(filter);
///     for (const BamRecord& record : query) {
///         // use record
///     }
/// \endcode
///
/// \note The whitelist versions only use Compare::EQUAL. Records will match at least one value
///       provided in the whitelist.
///
struct PbiZmwFilter : public internal::BasicDataFilterBase<int32_t, BasicLookupData::ZMW>
{
    PbiZmwFilter(const int32_t zmw, const Compare::Type cmp = Compare::EQUAL);
    PbiZmwFilter(const std::vector<int32_t>& whitelist);
    PbiZmwFilter(std::vector<int32_t>&& whitelist);
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/PbiFilterTypes.inl"

#endif // PBIFILTERTYPES_H
