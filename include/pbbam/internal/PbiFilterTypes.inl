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
/// \file PbiFilterTypes.inl
/// \brief Inline implementations for the built-in PBI filters.
//
// Author: Derek Barnett

#include "pbbam/PbiFilterTypes.h"

namespace PacBio {
namespace BAM {

namespace internal {

template <typename T>
inline FilterBase<T>::FilterBase(const T& value, const Compare::Type cmp)
    : value_(value)
    , cmp_(cmp)
{ }

template <typename T>
inline FilterBase<T>::FilterBase(T&& value, const Compare::Type cmp)
    : value_(std::move(value))
    , cmp_(cmp)
{ }

template <typename T>
inline FilterBase<T>::FilterBase(const std::vector<T>& values)
    : multiValue_(values)
{ }

template <typename T>
inline FilterBase<T>::FilterBase(std::vector<T>&& values)
    : multiValue_(std::move(values))
{ }

// BarcodeDataFilterBase

template<typename T, BarcodeLookupData::Field field>
inline BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(const T& value, const Compare::Type cmp)
    : FilterBase<T>(value, cmp)
{ }

template<typename T, BarcodeLookupData::Field field>
inline BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(T&& value, const Compare::Type cmp)
    : FilterBase<T>(std::move(value), cmp)
{ }

template<typename T, BarcodeLookupData::Field field>
inline BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(const std::vector<T>& values)
    : FilterBase<T>(values)
{ }

template<typename T, BarcodeLookupData::Field field>
inline BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(std::vector<T>&& values)
    : FilterBase<T>(std::move(values))
{ }

template<typename T, BarcodeLookupData::Field field>
inline IndexList BarcodeDataFilterBase<T, field>::Lookup(const PbiIndex& index) const
{
    if (FilterBase<T>::multiValue_ == boost::none)
        return index.BarcodeData().Indices(field, FilterBase<T>::value_, FilterBase<T>::cmp_);
    else
        return index.BarcodeData().IndicesMulti(field, FilterBase<T>::multiValue_.get());
}

// BasicDataFilterBase

template<typename T, BasicLookupData::Field field>
inline BasicDataFilterBase<T, field>::BasicDataFilterBase(const T& value, const Compare::Type cmp)
    : FilterBase<T>(value, cmp)
{ }

template<typename T, BasicLookupData::Field field>
inline BasicDataFilterBase<T, field>::BasicDataFilterBase(T&& value, const Compare::Type cmp)
    : FilterBase<T>(std::move(value), cmp)
{ }

template<typename T, BasicLookupData::Field field>
inline BasicDataFilterBase<T, field>::BasicDataFilterBase(const std::vector<T>& values)
    : FilterBase<T>(values)
{ }

template<typename T, BasicLookupData::Field field>
inline BasicDataFilterBase<T, field>::BasicDataFilterBase(std::vector<T>&& values)
    : FilterBase<T>(std::move(values))
{ }

template<typename T, BasicLookupData::Field field>
inline IndexList BasicDataFilterBase<T, field>::Lookup(const PbiIndex& index) const
{
    if (FilterBase<T>::multiValue_ == boost::none)
        return index.BasicData().Indices(field, FilterBase<T>::value_, FilterBase<T>::cmp_);
    else
        return index.BasicData().IndicesMulti(field, FilterBase<T>::multiValue_.get());
}

// MappedDataFilterBase

template<typename T, MappedLookupData::Field field>
inline MappedDataFilterBase<T, field>::MappedDataFilterBase(const T& value, const Compare::Type cmp)
    : FilterBase<T>(value, cmp)
{ }

template<typename T, MappedLookupData::Field field>
inline MappedDataFilterBase<T, field>::MappedDataFilterBase(T&& value, const Compare::Type cmp)
    : FilterBase<T>(std::move(value), cmp)
{ }

template<typename T, MappedLookupData::Field field>
inline MappedDataFilterBase<T, field>::MappedDataFilterBase(const std::vector<T>& values)
    : FilterBase<T>(values)
{ }

template<typename T, MappedLookupData::Field field>
inline MappedDataFilterBase<T, field>::MappedDataFilterBase(std::vector<T>&& values)
    : FilterBase<T>(std::move(values))
{ }

template<typename T, MappedLookupData::Field field>
inline IndexList MappedDataFilterBase<T, field>::Lookup(const PbiIndex& index) const
{
    if (FilterBase<T>::multiValue_ == boost::none)
        return index.MappedData().Indices(field, FilterBase<T>::value_, FilterBase<T>::cmp_);
    else
        return index.MappedData().IndicesMulti(field, FilterBase<T>::multiValue_.get());
}

} // namespace internal

// PbiAlignedEndFilter

inline PbiAlignedEndFilter::PbiAlignedEndFilter(const uint32_t position, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, MappedLookupData::A_END>(position, cmp)
{ }

// PbiAlignedLengthFilter

inline PbiAlignedLengthFilter::PbiAlignedLengthFilter(const uint32_t length, const Compare::Type cmp)
    : internal::FilterBase<uint32_t>(length, cmp)
{ }

// PbiAlignedStartFilter

inline PbiAlignedStartFilter::PbiAlignedStartFilter(const uint32_t position, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, MappedLookupData::A_START>(position, cmp)
{ }

// PbiAlignedStrandFilter

inline PbiAlignedStrandFilter::PbiAlignedStrandFilter(const Strand strand, const Compare::Type cmp)
    : internal::MappedDataFilterBase<Strand, MappedLookupData::STRAND>(strand, cmp)
{
    if (cmp != Compare::EQUAL && cmp != Compare::NOT_EQUAL) {
        auto msg = std::string{ "Compare type: " };
        msg += Compare::TypeToName(cmp);
        msg += " not supported for PbiAlignedStrandFilter (use one of Compare::EQUAL or Compare::NOT_EQUAL).";
        throw std::runtime_error(msg);
    }
}

// PbiBarcodeFilter

inline PbiBarcodeFilter::PbiBarcodeFilter(const uint16_t barcode, const Compare::Type cmp)
    : compositeFilter_{ PbiFilter::Union({ PbiBarcodeForwardFilter{barcode,cmp},
                                           PbiBarcodeReverseFilter{barcode,cmp}
                                         })
                      }
{ }

inline PbiBarcodeFilter::PbiBarcodeFilter(const std::vector<uint16_t> &whitelist)
    : compositeFilter_{ PbiFilter::Union({ PbiBarcodeForwardFilter{whitelist},
                                           PbiBarcodeReverseFilter{whitelist}
                                         })
                      }
{ }

inline PbiBarcodeFilter::PbiBarcodeFilter(std::vector<uint16_t> &&whitelist)
    : compositeFilter_{ PbiFilter::Union({ PbiBarcodeForwardFilter{std::move(whitelist)},
                                           PbiBarcodeReverseFilter{std::move(whitelist)}
                                         })
                      }
{ }

inline IndexList PbiBarcodeFilter::Lookup(const PbiIndex& idx) const
{ return compositeFilter_.Lookup(idx); }

// PbiBarcodeForwardFilter

inline PbiBarcodeForwardFilter::PbiBarcodeForwardFilter(const uint16_t bcFwdId, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_FORWARD>(bcFwdId, cmp)
{ }

inline PbiBarcodeForwardFilter::PbiBarcodeForwardFilter(const std::vector<uint16_t>& whitelist)
    : internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_FORWARD>(whitelist)
{ }

inline PbiBarcodeForwardFilter::PbiBarcodeForwardFilter(std::vector<uint16_t>&& whitelist)
    : internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_FORWARD>(std::move(whitelist))
{ }

// PbiBarcodeQualityFilter

inline PbiBarcodeQualityFilter::PbiBarcodeQualityFilter(const uint8_t bcQuality, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<uint8_t, BarcodeLookupData::BC_QUALITY>(bcQuality, cmp)
{ }

// PbiBarcodeReverseFilter

inline PbiBarcodeReverseFilter::PbiBarcodeReverseFilter(const uint16_t bcRevId, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_REVERSE>(bcRevId, cmp)
{ }

inline PbiBarcodeReverseFilter::PbiBarcodeReverseFilter(const std::vector<uint16_t>& whitelist)
    : internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_REVERSE>(whitelist)
{ }

inline PbiBarcodeReverseFilter::PbiBarcodeReverseFilter(std::vector<uint16_t>&& whitelist)
    : internal::BarcodeDataFilterBase<uint16_t, BarcodeLookupData::BC_REVERSE>(std::move(whitelist))
{ }

// PbiBarcodesFilter

inline PbiBarcodesFilter::PbiBarcodesFilter(const std::pair<uint16_t, uint16_t> barcodes, const Compare::Type cmp)
    : PbiBarcodesFilter(barcodes.first, barcodes.second, cmp)
{ }

inline PbiBarcodesFilter::PbiBarcodesFilter(const uint16_t bcForward, const uint16_t bcReverse, const Compare::Type cmp)
    : compositeFilter_{ PbiFilter::Intersection({ PbiBarcodeForwardFilter{bcForward,cmp},
                                                  PbiBarcodeReverseFilter{bcReverse,cmp}
                                                })
                      }
{ }

inline IndexList PbiBarcodesFilter::Lookup(const PbiIndex& idx) const
{ return compositeFilter_.Lookup(idx); }

// PbiIdentityFilter

inline PbiIdentityFilter::PbiIdentityFilter(const float identity,
                                            const Compare::Type cmp)
    : internal::FilterBase<float>(identity, cmp)
{ }

// PbiMapQualityFilter

inline PbiMapQualityFilter::PbiMapQualityFilter(const uint8_t mapQual, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint8_t, MappedLookupData::MAP_QUALITY>(mapQual, cmp)
{ }

// PbiMovieNameFilter

inline PbiMovieNameFilter::PbiMovieNameFilter(const std::string& movieName)
    : internal::FilterBase<std::string>(movieName, Compare::EQUAL)
{ }

inline PbiMovieNameFilter::PbiMovieNameFilter(const std::vector<std::string>& whitelist)
    : internal::FilterBase<std::string>(whitelist)
{ }

inline PbiMovieNameFilter::PbiMovieNameFilter(std::vector<std::string>&& whitelist)
    : internal::FilterBase<std::string>(std::move(whitelist))
{ }

// PbiNumDeletedBasesFilter

inline PbiNumDeletedBasesFilter::PbiNumDeletedBasesFilter(const size_t numDeletions, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, MappedLookupData::N_DEL>(numDeletions, cmp)
{ }

// PbiNumInsertedBasesFilter

inline PbiNumInsertedBasesFilter::PbiNumInsertedBasesFilter(const size_t numInsertions, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, MappedLookupData::N_INS>(numInsertions, cmp)
{ }

// PbiNumMatchesFilter

inline PbiNumMatchesFilter::PbiNumMatchesFilter(const size_t numMatchedBases, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, MappedLookupData::N_M>(numMatchedBases, cmp)
{ }

// PbiNumMismatchesFilter

inline PbiNumMismatchesFilter::PbiNumMismatchesFilter(const size_t numMismatchedBases, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, MappedLookupData::N_MM>(numMismatchedBases, cmp)
{ }

// PbiQueryEndFilter

inline PbiQueryEndFilter::PbiQueryEndFilter(const int32_t position, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::Q_END>(position, cmp)
{ }

// PbiQueryLengthFilter

inline PbiQueryLengthFilter::PbiQueryLengthFilter(const int32_t length, const Compare::Type cmp)
    : internal::FilterBase<int32_t>(length, cmp)
{ }

// PbiQueryNameFilter

inline PbiQueryNameFilter::PbiQueryNameFilter(const std::string& qname)
    : internal::FilterBase<std::string>(qname, Compare::EQUAL)
{ }

inline PbiQueryNameFilter::PbiQueryNameFilter(const std::vector<std::string>& whitelist)
    : internal::FilterBase<std::string>(whitelist)
{ }

inline PbiQueryNameFilter::PbiQueryNameFilter(std::vector<std::string>&& whitelist)
    : internal::FilterBase<std::string>(std::move(whitelist))
{ }

// PbiQueryStartFilter

inline PbiQueryStartFilter::PbiQueryStartFilter(const int32_t position, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::Q_START>(position, cmp)
{ }

// PbiReadAccuracyFilter

inline PbiReadAccuracyFilter::PbiReadAccuracyFilter(const Accuracy accuracy, const Compare::Type cmp)
    : internal::BasicDataFilterBase<Accuracy, BasicLookupData::READ_QUALITY>(accuracy, cmp)
{ }

// PbiReadGroupFilter

inline PbiReadGroupFilter::PbiReadGroupFilter(const int32_t rgId, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>(rgId, cmp)
{ }

inline PbiReadGroupFilter::PbiReadGroupFilter(const std::string rgId, const Compare::Type cmp)
    : PbiReadGroupFilter(ReadGroupInfo::IdToInt(rgId), cmp)
{ }

inline PbiReadGroupFilter::PbiReadGroupFilter(const ReadGroupInfo& rg, const Compare::Type cmp)
    : PbiReadGroupFilter(rg.Id(), cmp)
{ }

inline PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<int32_t>& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>(whitelist)
{ }

inline PbiReadGroupFilter::PbiReadGroupFilter(std::vector<int32_t>&& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>(std::move(whitelist))
{ }

inline PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<std::string>& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>(std::vector<int32_t>())
{
    multiValue_->reserve(whitelist.size());
    for (const auto& rg : whitelist)
        multiValue_->push_back(ReadGroupInfo::IdToInt(rg));
}

inline PbiReadGroupFilter::PbiReadGroupFilter(std::vector<std::string>&& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>(std::vector<int32_t>())
{
    multiValue_->reserve(whitelist.size());
    for (auto&& rg : whitelist)
        multiValue_->push_back(ReadGroupInfo::IdToInt(rg));
}

inline PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<ReadGroupInfo>& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>(std::vector<int32_t>())
{
    multiValue_->reserve(whitelist.size());
    for (const auto& rg : whitelist)
        multiValue_->push_back(ReadGroupInfo::IdToInt(rg.Id()));
}

inline PbiReadGroupFilter::PbiReadGroupFilter(std::vector<ReadGroupInfo>&& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::RG_ID>(std::vector<int32_t>())
{
    multiValue_->reserve(whitelist.size());
    for (auto&& rg : whitelist)
        multiValue_->push_back(ReadGroupInfo::IdToInt(rg.Id()));
}

// PbiReferenceEndFilter

inline PbiReferenceEndFilter::PbiReferenceEndFilter(const uint32_t tEnd, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, MappedLookupData::T_END>(tEnd, cmp)
{ }

// PbiReferenceIdFilter

inline PbiReferenceIdFilter::PbiReferenceIdFilter(const int32_t tId, const Compare::Type cmp)
    : internal::MappedDataFilterBase<int32_t, MappedLookupData::T_ID>(tId, cmp)
{ }

inline PbiReferenceIdFilter::PbiReferenceIdFilter(const std::vector<int32_t>& whitelist)
    : internal::MappedDataFilterBase<int32_t, MappedLookupData::T_ID>(whitelist)
{ }

inline PbiReferenceIdFilter::PbiReferenceIdFilter(std::vector<int32_t>&& whitelist)
    : internal::MappedDataFilterBase<int32_t, MappedLookupData::T_ID>(std::move(whitelist))
{ }

// PbiReferenceNameFilter

inline PbiReferenceNameFilter::PbiReferenceNameFilter(const std::string& rname,
                                                      const Compare::Type cmp)
    : internal::FilterBase<std::string>(rname, cmp)
{
    if (cmp != Compare::EQUAL && cmp != Compare::NOT_EQUAL) {
        auto msg = std::string{ "Compare type: " };
        msg += Compare::TypeToName(cmp);
        msg += " not supported for PbiReferenceNameFilter (use one of Compare::EQUAL or Compare::NOT_EQUAL).";
        throw std::runtime_error(msg);
    }
}

inline PbiReferenceNameFilter::PbiReferenceNameFilter(const std::vector<std::string>& whitelist)
    : internal::FilterBase<std::string>(whitelist)
{ }

inline PbiReferenceNameFilter::PbiReferenceNameFilter(std::vector<std::string>&& whitelist)
    : internal::FilterBase<std::string>(std::move(whitelist))
{ }

// PbiReferenceStartFilter

inline PbiReferenceStartFilter::PbiReferenceStartFilter(const uint32_t tStart, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, MappedLookupData::T_START>(tStart, cmp)
{ }

// PbiZmwFilter

inline PbiZmwFilter::PbiZmwFilter(const int32_t zmw, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::ZMW>(zmw, cmp)
{ }

inline PbiZmwFilter::PbiZmwFilter(const std::vector<int32_t>& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::ZMW>(whitelist)
{ }

inline PbiZmwFilter::PbiZmwFilter(std::vector<int32_t>&& whitelist)
    : internal::BasicDataFilterBase<int32_t, BasicLookupData::ZMW>(std::move(whitelist))
{ }

} // namespace BAM
} // namespace PacBio
