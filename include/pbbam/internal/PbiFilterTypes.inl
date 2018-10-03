// File Description
/// \file PbiFilterTypes.inl
/// \brief Inline implementations for the built-in PBI filters.
//
// Author: Derek Barnett

#include "pbbam/PbiFilterTypes.h"
#include <cassert>
#include <stdexcept>

#include <boost/functional/hash/hash.hpp>

namespace PacBio {
namespace BAM {

namespace internal {

template <typename T>
inline FilterBase<T>::FilterBase(T value, const Compare::Type cmp)
    : value_{std::move(value)}
    , cmp_{cmp}
{ }

template <typename T>
inline FilterBase<T>::FilterBase(std::vector<T> values, const Compare::Type cmp)
    : multiValue_{std::move(values)}
    , cmp_{cmp}
{ }

template<typename T>
inline bool FilterBase<T>::CompareHelper(const T& lhs) const
{
    if (multiValue_ == boost::none)
        return CompareSingleHelper(lhs);
    else
        return CompareMultiHelper(lhs);
}

template<typename T>
inline bool FilterBase<T>::CompareMultiHelper(const T& lhs) const
{
    if (cmp_ == Compare::EQUAL)
    {
        // check provided value against all filter criteria,
        // return true on any exact match
        auto iter = multiValue_.get().cbegin();
        const auto end  = multiValue_.get().cend();
        for (; iter != end; ++iter) {
            if (*iter == lhs)
                return true;
        }
        return false; // no matches
    }
    else if (cmp_ == Compare::NOT_EQUAL)
    {
        // check provided value against all filter criteria,
        // return true on any exact match
        auto iter = multiValue_.get().cbegin();
        const auto end  = multiValue_.get().cend();
        for (; iter != end; ++iter) {
            if (*iter == lhs)
                return false;
        }
        return true;
    }
    else
        throw std::runtime_error{"unsupported compare type on multivalue filter"};
}

template<typename T>
inline bool FilterBase<T>::CompareSingleHelper(const T& lhs) const
{
    return Compare::Check(lhs, value_, cmp_);
}

template<>
inline bool FilterBase<LocalContextFlags>::CompareSingleHelper(const LocalContextFlags& lhs) const
{
    switch(cmp_) {
        case Compare::EQUAL:              return lhs == value_;
        case Compare::LESS_THAN:          return lhs < value_;
        case Compare::LESS_THAN_EQUAL:    return lhs <= value_;
        case Compare::GREATER_THAN:       return lhs > value_;
        case Compare::GREATER_THAN_EQUAL: return lhs >= value_;
        case Compare::NOT_EQUAL:          return lhs != value_;
        case Compare::CONTAINS:           return ((lhs & value_) != 0);
        case Compare::NOT_CONTAINS:       return ((lhs & value_) == 0);

        default:
            assert(false);
            throw std::runtime_error{"unsupported compare type requested"};
    }
}

// BarcodeDataFilterBase

template<typename T, PbiFile::BarcodeField field>
inline BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(T value, const Compare::Type cmp)
    : FilterBase<T>{std::move(value), cmp}
{ }

template<typename T, PbiFile::BarcodeField field>
inline BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(std::vector<T> values, const Compare::Type cmp)
    : FilterBase<T>{std::move(values), cmp}
{ }

template<typename T, PbiFile::BarcodeField field>
inline bool BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase::Accepts(const PbiRawData& idx,
                                           const size_t row) const
{
    const PbiRawBarcodeData& barcodeData = idx.BarcodeData();
    switch (field) {
        case PbiFile::BarcodeField::BC_FORWARD: return FilterBase<T>::CompareHelper(barcodeData.bcForward_.at(row));
        case PbiFile::BarcodeField::BC_REVERSE: return FilterBase<T>::CompareHelper(barcodeData.bcReverse_.at(row));
        case PbiFile::BarcodeField::BC_QUALITY: return FilterBase<T>::CompareHelper(barcodeData.bcQual_.at(row));
        default:
            assert(false);
            throw std::runtime_error{"unsupported BarcodeData field requested"};
    }
}

// BasicDataFilterBase

template<typename T, PbiFile::BasicField field>
inline BasicDataFilterBase<T, field>::BasicDataFilterBase(T value, const Compare::Type cmp)
    : FilterBase<T>{std::move(value), cmp}
{ }

template<typename T, PbiFile::BasicField field>
inline BasicDataFilterBase<T, field>::BasicDataFilterBase(std::vector<T> values, const Compare::Type cmp)
    : FilterBase<T>{std::move(values), cmp}
{ }

template<typename T, PbiFile::BasicField field>
inline bool BasicDataFilterBase<T, field>::BasicDataFilterBase::Accepts(const PbiRawData& idx,
                                                                        const size_t row) const
{
    const PbiRawBasicData& basicData = idx.BasicData();
    switch (field) {
        case PbiFile::BasicField::RG_ID:        return FilterBase<T>::CompareHelper(basicData.rgId_.at(row));
        case PbiFile::BasicField::Q_START:      return FilterBase<T>::CompareHelper(basicData.qStart_.at(row));
        case PbiFile::BasicField::Q_END:        return FilterBase<T>::CompareHelper(basicData.qEnd_.at(row));
        case PbiFile::BasicField::ZMW:          return FilterBase<T>::CompareHelper(basicData.holeNumber_.at(row));
        case PbiFile::BasicField::READ_QUALITY: return FilterBase<T>::CompareHelper(basicData.readQual_.at(row));
        //   PbiFile::BasicField::CONTEXT_FLAG has its own specialization
        default:
            assert(false);
            throw std::runtime_error{"unsupported BasicData field requested"};
    }
}

// this typedef exists purely so that the next method signature isn't 2 screen widths long
using LocalContextFilter__ = BasicDataFilterBase<LocalContextFlags, PbiFile::BasicField::CONTEXT_FLAG>;

template<>
inline bool LocalContextFilter__::BasicDataFilterBase::Accepts(const PbiRawData& idx,
                                                               const size_t row) const
{
    const auto& basicData = idx.BasicData();
    const auto rowFlags = static_cast<LocalContextFlags>(basicData.ctxtFlag_.at(row));
    return FilterBase<LocalContextFlags>::CompareHelper(rowFlags);
}

// BasicDataFilterBase

template<typename T, PbiFile::MappedField field>
inline MappedDataFilterBase<T, field>::MappedDataFilterBase(T value, const Compare::Type cmp)
    : FilterBase<T>{std::move(value), cmp}
{ }

template<typename T, PbiFile::MappedField field>
inline MappedDataFilterBase<T, field>::MappedDataFilterBase(std::vector<T> values, const Compare::Type cmp)
    : FilterBase<T>{std::move(values), cmp}
{ }

template<>
inline bool MappedDataFilterBase<Strand, PbiFile::MappedField::STRAND>::MappedDataFilterBase::Accepts(const PbiRawData& idx,
                                                                                                  const size_t row) const
{
    const PbiRawMappedData& mappedData = idx.MappedData();
    const Strand strand = (mappedData.revStrand_.at(row) == 1 ? Strand::REVERSE : Strand::FORWARD);
    return FilterBase<Strand>::CompareHelper(strand);
}

template<typename T, PbiFile::MappedField field>
inline bool MappedDataFilterBase<T, field>::MappedDataFilterBase::Accepts(const PbiRawData& idx,
                                                                          const size_t row) const
{
    const PbiRawMappedData& mappedData = idx.MappedData();
    switch (field) {
        case PbiFile::MappedField::T_ID:        return FilterBase<T>::CompareHelper(mappedData.tId_.at(row));
        case PbiFile::MappedField::T_START:     return FilterBase<T>::CompareHelper(mappedData.tStart_.at(row));
        case PbiFile::MappedField::T_END:       return FilterBase<T>::CompareHelper(mappedData.tEnd_.at(row));
        case PbiFile::MappedField::A_START:     return FilterBase<T>::CompareHelper(mappedData.aStart_.at(row));
        case PbiFile::MappedField::A_END:       return FilterBase<T>::CompareHelper(mappedData.aEnd_.at(row));
        case PbiFile::MappedField::N_M:         return FilterBase<T>::CompareHelper(mappedData.nM_.at(row));
        case PbiFile::MappedField::N_MM:        return FilterBase<T>::CompareHelper(mappedData.nMM_.at(row));
        case PbiFile::MappedField::N_DEL:       return FilterBase<T>::CompareHelper(mappedData.NumDeletedBasesAt(row));
        case PbiFile::MappedField::N_INS:       return FilterBase<T>::CompareHelper(mappedData.NumInsertedBasesAt(row));
        case PbiFile::MappedField::MAP_QUALITY: return FilterBase<T>::CompareHelper(mappedData.mapQV_.at(row));
        default:
            assert(false);
            throw std::runtime_error{"unsupported MappedData field requested"};
    }
}

} // namespace internal

// PbiAlignedEndFilter

inline PbiAlignedEndFilter::PbiAlignedEndFilter(const uint32_t position, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::A_END>{position, cmp}
{ }

// PbiAlignedLengthFilter

inline PbiAlignedLengthFilter::PbiAlignedLengthFilter(const uint32_t length, const Compare::Type cmp)
    : internal::FilterBase<uint32_t>{length, cmp}
{ }

// PbiAlignedStartFilter

inline PbiAlignedStartFilter::PbiAlignedStartFilter(const uint32_t position, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::A_START>{position, cmp}
{ }

// PbiAlignedStrandFilter

inline PbiAlignedStrandFilter::PbiAlignedStrandFilter(const Strand strand, const Compare::Type cmp)
    : internal::MappedDataFilterBase<Strand, PbiFile::MappedField::STRAND>{strand, cmp}
{
    if (cmp != Compare::EQUAL && cmp != Compare::NOT_EQUAL) {
        throw std::runtime_error{"Compare type: " + Compare::TypeToName(cmp) + " not supported for PbiAlignedStrandFilter (use one of Compare::EQUAL or Compare::NOT_EQUAL)."};
    }
}

// PbiBarcodeFilter

inline PbiBarcodeFilter::PbiBarcodeFilter(const int16_t barcode, const Compare::Type cmp)
    : compositeFilter_{ PbiFilter::Union({ PbiBarcodeForwardFilter{barcode,cmp},
                                           PbiBarcodeReverseFilter{barcode,cmp}
                                         })
                      }
{ }

inline PbiBarcodeFilter::PbiBarcodeFilter(std::vector<int16_t> whitelist, const Compare::Type cmp)
    : compositeFilter_{ PbiFilter::Union({ PbiBarcodeForwardFilter{std::move(whitelist), cmp},
                                           PbiBarcodeReverseFilter{std::move(whitelist), cmp}
                                         })
                      }
{ }

inline bool PbiBarcodeFilter::Accepts(const PbiRawData& idx, const size_t row) const
{ return compositeFilter_.Accepts(idx, row); }

// PbiBarcodeForwardFilter

inline PbiBarcodeForwardFilter::PbiBarcodeForwardFilter(const int16_t bcFwdId, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_FORWARD>{bcFwdId, cmp}
{ }

inline PbiBarcodeForwardFilter::PbiBarcodeForwardFilter(std::vector<int16_t> whitelist, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_FORWARD>{std::move(whitelist), cmp}
{ }

// PbiBarcodeQualityFilter

inline PbiBarcodeQualityFilter::PbiBarcodeQualityFilter(const uint8_t bcQuality, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<uint8_t, PbiFile::BarcodeField::BC_QUALITY>{bcQuality, cmp}
{ }

// PbiBarcodeReverseFilter

inline PbiBarcodeReverseFilter::PbiBarcodeReverseFilter(const int16_t bcRevId, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_REVERSE>{bcRevId, cmp}
{ }

inline PbiBarcodeReverseFilter::PbiBarcodeReverseFilter(std::vector<int16_t> whitelist, const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_REVERSE>{std::move(whitelist), cmp}
{ }

// PbiBarcodesFilter

inline PbiBarcodesFilter::PbiBarcodesFilter(const std::pair<int16_t, int16_t> barcodes, const Compare::Type cmp)
    : PbiBarcodesFilter{barcodes.first, barcodes.second, cmp}
{ }

inline PbiBarcodesFilter::PbiBarcodesFilter(const int16_t bcForward, const int16_t bcReverse, const Compare::Type cmp)
    : compositeFilter_{ PbiFilter::Intersection({ PbiBarcodeForwardFilter{bcForward,cmp},
                                                  PbiBarcodeReverseFilter{bcReverse,cmp}
                                                })
                      }
{ }

inline bool PbiBarcodesFilter::Accepts(const PbiRawData& idx, const size_t row) const
{ return compositeFilter_.Accepts(idx, row); }

// PbiIdentityFilter

inline PbiIdentityFilter::PbiIdentityFilter(const float identity,
                                            const Compare::Type cmp)
    : internal::FilterBase<float>{identity, cmp}
{ }

// PbiLocalContextFilter

inline PbiLocalContextFilter::PbiLocalContextFilter(const LocalContextFlags& flags,
                                                    const Compare::Type cmp)
    : internal::BasicDataFilterBase<LocalContextFlags, PbiFile::BasicField::CONTEXT_FLAG>{flags, cmp}
{ }

// PbiMapQualityFilter

inline PbiMapQualityFilter::PbiMapQualityFilter(const uint8_t mapQual, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint8_t, PbiFile::MappedField::MAP_QUALITY>{mapQual, cmp}
{ }

// PbiMovieNameFilter

inline bool PbiMovieNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const bool found = compositeFilter_.Accepts(idx, row);
    if (cmp_ == Compare::EQUAL) return found;
    else if (cmp_ == Compare::NOT_EQUAL) return !found;
    else throw std::runtime_error{"unsupported compare type on movie name filter"};
}

// PbiNumDeletedBasesFilter

inline PbiNumDeletedBasesFilter::PbiNumDeletedBasesFilter(const size_t numDeletions, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_DEL>{numDeletions, cmp}
{ }

// PbiNumInsertedBasesFilter

inline PbiNumInsertedBasesFilter::PbiNumInsertedBasesFilter(const size_t numInsertions, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_INS>{numInsertions, cmp}
{ }

// PbiNumMatchesFilter

inline PbiNumMatchesFilter::PbiNumMatchesFilter(const size_t numMatchedBases, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_M>{numMatchedBases, cmp}
{ }

// PbiNumMismatchesFilter

inline PbiNumMismatchesFilter::PbiNumMismatchesFilter(const size_t numMismatchedBases, const Compare::Type cmp)
    : internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_MM>{numMismatchedBases, cmp}
{ }

// PbiQueryEndFilter

inline PbiQueryEndFilter::PbiQueryEndFilter(const int32_t position, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::Q_END>{position, cmp}
{ }

// PbiQueryLengthFilter

inline PbiQueryLengthFilter::PbiQueryLengthFilter(const int32_t length, const Compare::Type cmp)
    : internal::FilterBase<int32_t>{length, cmp}
{ }

// PbiQueryStartFilter

inline PbiQueryStartFilter::PbiQueryStartFilter(const int32_t position, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::Q_START>{position, cmp}
{ }

// PbiReadAccuracyFilter

inline PbiReadAccuracyFilter::PbiReadAccuracyFilter(const Accuracy accuracy, const Compare::Type cmp)
    : internal::BasicDataFilterBase<Accuracy, PbiFile::BasicField::READ_QUALITY>{accuracy, cmp}
{ }

// PbiReferenceEndFilter

inline PbiReferenceEndFilter::PbiReferenceEndFilter(const uint32_t tEnd, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::T_END>{tEnd, cmp}
{ }

// PbiReferenceIdFilter

inline PbiReferenceIdFilter::PbiReferenceIdFilter(const int32_t tId, const Compare::Type cmp)
    : internal::MappedDataFilterBase<int32_t, PbiFile::MappedField::T_ID>{tId, cmp}
{ }

inline PbiReferenceIdFilter::PbiReferenceIdFilter(std::vector<int32_t> whitelist, const Compare::Type cmp)
    : internal::MappedDataFilterBase<int32_t, PbiFile::MappedField::T_ID>{std::move(whitelist), cmp}
{ }

// PbiReferenceStartFilter

inline PbiReferenceStartFilter::PbiReferenceStartFilter(const uint32_t tStart, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::T_START>{tStart, cmp}
{ }

// PbiZmwFilter

inline PbiZmwFilter::PbiZmwFilter(const int32_t zmw, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::ZMW>{zmw, cmp}
{ }

inline PbiZmwFilter::PbiZmwFilter(std::vector<int32_t> whitelist, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::ZMW>{std::move(whitelist), cmp}
{ }

// PbiZmwModuloFilter

inline PbiZmwModuloFilter::PbiZmwModuloFilter(
        const uint32_t denominator,
        const uint32_t value,
        const FilterHash hashType,
        const Compare::Type cmp)
    : denominator_{denominator}
    , value_{value}
    , hash_{hashType}
    , cmp_{cmp}
{ }

inline uint32_t UnsignedLongIntCast(const int32_t zm)
{
    return static_cast<uint32_t>(zm);
}

inline uint32_t BoostHashCombine(const int32_t zm)
{
    constexpr static const uint16_t mask = 0xFFFF;

    const uint16_t upper = (zm >> 16) & mask;
    const uint16_t lower = zm & mask;

    // FIXME: discrepancies with Python API. Will return to nail down.

    size_t seed = 0;
    boost::hash_combine(seed, upper);
    boost::hash_combine(seed, lower);
    return static_cast<uint32_t>(seed);
}

inline bool PbiZmwModuloFilter::Accepts(const PbiRawData& idx,
                                        const size_t row) const
{
    const auto zm = idx.BasicData().holeNumber_.at(row);

    uint32_t hashValue;
    switch(hash_)
    {
        case FilterHash::UNSIGNED_LONG_CAST :
        {
            hashValue = UnsignedLongIntCast(zm);
            break;
        }

        case FilterHash::BOOST_HASH_COMBINE :
        {
            hashValue = BoostHashCombine(zm);
            break;
        }

        default:
            throw std::runtime_error{"unsupported filter hash type"};
    }

    const auto modResult = hashValue % denominator_;
    return Compare::Check(modResult, value_, cmp_);
}

} // namespace BAM
} // namespace PacBio
