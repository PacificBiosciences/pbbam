#ifndef PBBAM_PBIFILTERTYPES_INL
#define PBBAM_PBIFILTERTYPES_INL

#include <pbbam/Config.h>

#include <pbbam/PbiFilterTypes.h>

#include <boost/functional/hash/hash.hpp>

#include <stdexcept>

#include <cassert>

namespace PacBio {
namespace BAM {

namespace internal {

template <typename T>
FilterBase<T>::FilterBase(T value, const Compare::Type cmp) : value_{std::move(value)}, cmp_{cmp}
{
}

template <typename T>
FilterBase<T>::FilterBase(std::vector<T> values, const Compare::Type cmp)
    : multiValue_{std::move(values)}, cmp_{cmp}
{
    // "=="/"!=" can come in from XML, e.g. <Property Name="zmw" Operator="==" Value="(x,y,z)" />"
    // switch to whitelist/blacklist containment for multi-value filters
    if (cmp_ == Compare::EQUAL) {
        cmp_ = Compare::CONTAINS;
    } else if (cmp_ == Compare::NOT_EQUAL) {
        cmp_ = Compare::NOT_CONTAINS;
    }

    if (cmp_ != Compare::CONTAINS && cmp_ != Compare::NOT_CONTAINS) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: multi-valued filters (e.g. whitelists) can only check for "
            "containment."};
    }
}

template <typename T>
bool FilterBase<T>::CompareHelper(const T& lhs) const
{
    if (multiValue_) {
        return CompareMultiHelper(lhs);
    } else {
        return CompareSingleHelper(lhs);
    }
}

template <typename T>
bool FilterBase<T>::CompareMultiHelper(const T& lhs) const
{
    // multi-value filters are whitelist/blacklist
    assert(cmp_ == Compare::CONTAINS || cmp_ == Compare::NOT_CONTAINS);

    // whitelist - return true on any hit
    if (cmp_ == Compare::CONTAINS) {
        for (const auto& x : *multiValue_) {
            if (x == lhs) {
                return true;
            }
        }
        return false;
    }
    // blacklist - return false on any hit
    else {
        for (const auto& x : *multiValue_) {
            if (x == lhs) {
                return false;
            }
        }
        return true;
    }
}

template <typename T>
bool FilterBase<T>::CompareSingleHelper(const T& lhs) const
{
    return Compare::Check(lhs, value_, cmp_);
}

template <>
inline bool FilterBase<Data::LocalContextFlags>::CompareSingleHelper(const Data::LocalContextFlags& lhs) const
{
    switch (cmp_) {
        case Compare::EQUAL:
            return lhs == value_;
        case Compare::LESS_THAN:
            return lhs < value_;
        case Compare::LESS_THAN_EQUAL:
            return lhs <= value_;
        case Compare::GREATER_THAN:
            return lhs > value_;
        case Compare::GREATER_THAN_EQUAL:
            return lhs >= value_;
        case Compare::NOT_EQUAL:
            return lhs != value_;
        case Compare::CONTAINS:
            return ((lhs & value_) != 0);
        case Compare::NOT_CONTAINS:
            return ((lhs & value_) == 0);

        default:
            assert(false);
            throw std::runtime_error{"[pbbam] PBI filter ERROR: unknown compare type '" +
                                     Compare::TypeToName(cmp_) + "'"};
    }
}

// BarcodeDataFilterBase

template <typename T, PbiFile::BarcodeField field>
BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(T value, const Compare::Type cmp)
    : FilterBase<T>{std::move(value), cmp}
{
}

template <typename T, PbiFile::BarcodeField field>
BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase(std::vector<T> values,
                                                       const Compare::Type cmp)
    : FilterBase<T>{std::move(values), cmp}
{
}

template <typename T, PbiFile::BarcodeField field>
bool BarcodeDataFilterBase<T, field>::BarcodeDataFilterBase::Accepts(const PbiRawData& idx,
                                                                     const std::size_t row) const
{
    const PbiRawBarcodeData& barcodeData = idx.BarcodeData();
    switch (field) {
        case PbiFile::BarcodeField::BC_FORWARD:
            return FilterBase<T>::CompareHelper(barcodeData.bcForward_.at(row));
        case PbiFile::BarcodeField::BC_REVERSE:
            return FilterBase<T>::CompareHelper(barcodeData.bcReverse_.at(row));
        case PbiFile::BarcodeField::BC_QUALITY:
            return FilterBase<T>::CompareHelper(barcodeData.bcQual_.at(row));
        default:
            assert(false);
            throw std::runtime_error{"[pbbam] PBI filter ERROR: unknown barcode field requested."};
    }
}

// BasicDataFilterBase

template <typename T, PbiFile::BasicField field>
BasicDataFilterBase<T, field>::BasicDataFilterBase(T value, const Compare::Type cmp)
    : FilterBase<T>{std::move(value), cmp}
{
}

template <typename T, PbiFile::BasicField field>
BasicDataFilterBase<T, field>::BasicDataFilterBase(std::vector<T> values, const Compare::Type cmp)
    : FilterBase<T>{std::move(values), cmp}
{
}

template <typename T, PbiFile::BasicField field>
bool BasicDataFilterBase<T, field>::BasicDataFilterBase::Accepts(const PbiRawData& idx,
                                                                 const std::size_t row) const
{
    const PbiRawBasicData& basicData = idx.BasicData();
    switch (field) {
        case PbiFile::BasicField::RG_ID:
            return FilterBase<T>::CompareHelper(basicData.rgId_.at(row));
        case PbiFile::BasicField::Q_START:
            return FilterBase<T>::CompareHelper(basicData.qStart_.at(row));
        case PbiFile::BasicField::Q_END:
            return FilterBase<T>::CompareHelper(basicData.qEnd_.at(row));
        case PbiFile::BasicField::ZMW:
            return FilterBase<T>::CompareHelper(basicData.holeNumber_.at(row));
        case PbiFile::BasicField::READ_QUALITY:
            return FilterBase<T>::CompareHelper(basicData.readQual_.at(row));
        // NOTE(DB): PbiFile::BasicField::CONTEXT_FLAG has its own specialization
        default:
            assert(false);
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unknown basic data field requested."};
    }
}

// this typedef exists purely so that the next method signature isn't 2 screen widths long
using LocalContextFilterInternal =
    BasicDataFilterBase<Data::LocalContextFlags, PbiFile::BasicField::CONTEXT_FLAG>;

template <>
inline bool LocalContextFilterInternal::BasicDataFilterBase::Accepts(const PbiRawData& idx,
                                                               const std::size_t row) const
{
    const auto& basicData = idx.BasicData();
    const auto rowFlags = static_cast<Data::LocalContextFlags>(basicData.ctxtFlag_.at(row));
    return FilterBase<Data::LocalContextFlags>::CompareHelper(rowFlags);
}

// BasicDataFilterBase

template <typename T, PbiFile::MappedField field>
MappedDataFilterBase<T, field>::MappedDataFilterBase(T value, const Compare::Type cmp)
    : FilterBase<T>{std::move(value), cmp}
{
}

template <typename T, PbiFile::MappedField field>
MappedDataFilterBase<T, field>::MappedDataFilterBase(std::vector<T> values, const Compare::Type cmp)
    : FilterBase<T>{std::move(values), cmp}
{
}

template <>
inline bool
MappedDataFilterBase<Data::Strand, PbiFile::MappedField::STRAND>::MappedDataFilterBase::Accepts(
    const PbiRawData& idx, const std::size_t row) const
{
    const PbiRawMappedData& mappedData = idx.MappedData();
    const Data::Strand strand = (mappedData.revStrand_.at(row) == 1 ? Data::Strand::REVERSE : Data::Strand::FORWARD);
    return FilterBase<Data::Strand>::CompareHelper(strand);
}

template <typename T, PbiFile::MappedField field>
bool MappedDataFilterBase<T, field>::MappedDataFilterBase::Accepts(const PbiRawData& idx,
                                                                   const std::size_t row) const
{
    const PbiRawMappedData& mappedData = idx.MappedData();
    switch (field) {
        case PbiFile::MappedField::T_ID:
            return FilterBase<T>::CompareHelper(mappedData.tId_.at(row));
        case PbiFile::MappedField::T_START:
            return FilterBase<T>::CompareHelper(mappedData.tStart_.at(row));
        case PbiFile::MappedField::T_END:
            return FilterBase<T>::CompareHelper(mappedData.tEnd_.at(row));
        case PbiFile::MappedField::A_START:
            return FilterBase<T>::CompareHelper(mappedData.aStart_.at(row));
        case PbiFile::MappedField::A_END:
            return FilterBase<T>::CompareHelper(mappedData.aEnd_.at(row));
        case PbiFile::MappedField::N_M:
            return FilterBase<T>::CompareHelper(mappedData.nM_.at(row));
        case PbiFile::MappedField::N_MM:
            return FilterBase<T>::CompareHelper(mappedData.nMM_.at(row));
        case PbiFile::MappedField::N_DEL:
            return FilterBase<T>::CompareHelper(mappedData.NumDeletedBasesAt(row));
        case PbiFile::MappedField::N_INS:
            return FilterBase<T>::CompareHelper(mappedData.NumInsertedBasesAt(row));
        case PbiFile::MappedField::MAP_QUALITY:
            return FilterBase<T>::CompareHelper(mappedData.mapQV_.at(row));
        default:
            assert(false);
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unknown mapped data field requested."};
    }
}

}  // namespace internal

// PbiAlignedEndFilter

inline PbiAlignedEndFilter::PbiAlignedEndFilter(const uint32_t position, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::A_END>{position, cmp}
{
}

// PbiAlignedLengthFilter

inline PbiAlignedLengthFilter::PbiAlignedLengthFilter(const uint32_t length,
                                                      const Compare::Type cmp)
    : internal::FilterBase<uint32_t>{length, cmp}
{
}

// PbiAlignedStartFilter

inline PbiAlignedStartFilter::PbiAlignedStartFilter(const uint32_t position,
                                                    const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::A_START>{position, cmp}
{
}

// PbiAlignedStrandFilter

inline PbiAlignedStrandFilter::PbiAlignedStrandFilter(const Data::Strand strand, const Compare::Type cmp)
    : internal::MappedDataFilterBase<Data::Strand, PbiFile::MappedField::STRAND>{strand, cmp}
{
    if (cmp != Compare::EQUAL && cmp != Compare::NOT_EQUAL) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: compare type for aligned strand must be either EQUAL or "
            "NOT_EQUAL"};
    }
}

// PbiBarcodeFilter

inline PbiBarcodeFilter::PbiBarcodeFilter(const int16_t barcode, const Compare::Type cmp)
    : compositeFilter_{PbiFilter::Union(
          {PbiBarcodeForwardFilter{barcode, cmp}, PbiBarcodeReverseFilter{barcode, cmp}})}
{
}

inline PbiBarcodeFilter::PbiBarcodeFilter(std::vector<int16_t> barcodes, const Compare::Type cmp)
    : compositeFilter_{PbiFilter::Union(
          {PbiBarcodeForwardFilter{barcodes, cmp}, PbiBarcodeReverseFilter{barcodes, cmp}})}
{
}

inline bool PbiBarcodeFilter::Accepts(const PbiRawData& idx, const std::size_t row) const
{
    return compositeFilter_.Accepts(idx, row);
}

// PbiBarcodeForwardFilter

inline PbiBarcodeForwardFilter::PbiBarcodeForwardFilter(const int16_t bcFwdId,
                                                        const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_FORWARD>{bcFwdId, cmp}
{
}

inline PbiBarcodeForwardFilter::PbiBarcodeForwardFilter(std::vector<int16_t> barcodes,
                                                        const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_FORWARD>{
          std::move(barcodes), cmp}
{
}

// PbiBarcodeQualityFilter

inline PbiBarcodeQualityFilter::PbiBarcodeQualityFilter(const uint8_t bcQuality,
                                                        const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<uint8_t, PbiFile::BarcodeField::BC_QUALITY>{bcQuality, cmp}
{
}

// PbiBarcodeReverseFilter

inline PbiBarcodeReverseFilter::PbiBarcodeReverseFilter(const int16_t bcRevId,
                                                        const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_REVERSE>{bcRevId, cmp}
{
}

inline PbiBarcodeReverseFilter::PbiBarcodeReverseFilter(std::vector<int16_t> barcodes,
                                                        const Compare::Type cmp)
    : internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_REVERSE>{
          std::move(barcodes), cmp}
{
}

// PbiBarcodesFilter

inline PbiBarcodesFilter::PbiBarcodesFilter(const std::pair<int16_t, int16_t> barcodes,
                                            const Compare::Type cmp)
    : PbiBarcodesFilter{barcodes.first, barcodes.second, cmp}
{
}

inline PbiBarcodesFilter::PbiBarcodesFilter(const int16_t bcForward, const int16_t bcReverse,
                                            const Compare::Type cmp)
    : compositeFilter_{PbiFilter::Intersection(
          {PbiBarcodeForwardFilter{bcForward, cmp}, PbiBarcodeReverseFilter{bcReverse, cmp}})}
{
}

inline bool PbiBarcodesFilter::Accepts(const PbiRawData& idx, const std::size_t row) const
{
    return compositeFilter_.Accepts(idx, row);
}

// PbiIdentityFilter

inline PbiIdentityFilter::PbiIdentityFilter(const float identity, const Compare::Type cmp)
    : internal::FilterBase<float>{identity, cmp}
{
}

// PbiLocalContextFilter

inline PbiLocalContextFilter::PbiLocalContextFilter(const Data::LocalContextFlags& flags,
                                                    const Compare::Type cmp)
    : internal::BasicDataFilterBase<Data::LocalContextFlags, PbiFile::BasicField::CONTEXT_FLAG>{flags,
                                                                                          cmp}
{
}

// PbiMapQualityFilter

inline PbiMapQualityFilter::PbiMapQualityFilter(const uint8_t mapQual, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint8_t, PbiFile::MappedField::MAP_QUALITY>{mapQual, cmp}
{
}

// PbiNumDeletedBasesFilter

inline PbiNumDeletedBasesFilter::PbiNumDeletedBasesFilter(const std::size_t numDeletions,
                                                          const Compare::Type cmp)
    : internal::MappedDataFilterBase<std::size_t, PbiFile::MappedField::N_DEL>{numDeletions, cmp}
{
}

// PbiNumInsertedBasesFilter

inline PbiNumInsertedBasesFilter::PbiNumInsertedBasesFilter(const std::size_t numInsertions,
                                                            const Compare::Type cmp)
    : internal::MappedDataFilterBase<std::size_t, PbiFile::MappedField::N_INS>{numInsertions, cmp}
{
}

// PbiNumMatchesFilter

inline PbiNumMatchesFilter::PbiNumMatchesFilter(const std::size_t numMatchedBases,
                                                const Compare::Type cmp)
    : internal::MappedDataFilterBase<std::size_t, PbiFile::MappedField::N_M>{numMatchedBases, cmp}
{
}

// PbiNumMismatchesFilter

inline PbiNumMismatchesFilter::PbiNumMismatchesFilter(const std::size_t numMismatchedBases,
                                                      const Compare::Type cmp)
    : internal::MappedDataFilterBase<std::size_t, PbiFile::MappedField::N_MM>{numMismatchedBases, cmp}
{
}

// PbiQueryEndFilter

inline PbiQueryEndFilter::PbiQueryEndFilter(const int32_t position, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::Q_END>{position, cmp}
{
}

// PbiQueryLengthFilter

inline PbiQueryLengthFilter::PbiQueryLengthFilter(const int32_t length, const Compare::Type cmp)
    : internal::FilterBase<int32_t>{length, cmp}
{
}

// PbiQueryStartFilter

inline PbiQueryStartFilter::PbiQueryStartFilter(const int32_t position, const Compare::Type cmp)
    : internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::Q_START>{position, cmp}
{
}

// PbiReadAccuracyFilter

inline PbiReadAccuracyFilter::PbiReadAccuracyFilter(const Data::Accuracy accuracy,
                                                    const Compare::Type cmp)
    : internal::BasicDataFilterBase<Data::Accuracy, PbiFile::BasicField::READ_QUALITY>{accuracy, cmp}
{
}

// PbiReferenceEndFilter

inline PbiReferenceEndFilter::PbiReferenceEndFilter(const uint32_t tEnd, const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::T_END>{tEnd, cmp}
{
}

// PbiReferenceIdFilter

inline PbiReferenceIdFilter::PbiReferenceIdFilter(const int32_t tId, const Compare::Type cmp)
    : internal::MappedDataFilterBase<int32_t, PbiFile::MappedField::T_ID>{tId, cmp}
{
}

inline PbiReferenceIdFilter::PbiReferenceIdFilter(std::vector<int32_t> tIds,
                                                  const Compare::Type cmp)
    : internal::MappedDataFilterBase<int32_t, PbiFile::MappedField::T_ID>{std::move(tIds), cmp}
{
}

// PbiReferenceStartFilter

inline PbiReferenceStartFilter::PbiReferenceStartFilter(const uint32_t tStart,
                                                        const Compare::Type cmp)
    : internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::T_START>{tStart, cmp}
{
}

// PbiZmwModuloFilter

inline PbiZmwModuloFilter::PbiZmwModuloFilter(const uint32_t denominator, const uint32_t value,
                                              const FilterHash hashType, const Compare::Type cmp)
    : denominator_{denominator}, value_{value}, hash_{hashType}, cmp_{cmp}
{
}

inline uint32_t UnsignedLongIntCast(const int32_t zm) { return static_cast<uint32_t>(zm); }

inline uint32_t BoostHashCombine(const int32_t zm)
{
    constexpr uint16_t MASK = 0xFFFF;

    const uint16_t upper = (zm >> 16) & MASK;
    const uint16_t lower = zm & MASK;

    // FIXME: discrepancies with Python API. Will return to nail down.

    std::size_t seed = 0;
    boost::hash_combine(seed, upper);
    boost::hash_combine(seed, lower);
    return static_cast<uint32_t>(seed);
}

inline bool PbiZmwModuloFilter::Accepts(const PbiRawData& idx, const std::size_t row) const
{
    const auto zm = idx.BasicData().holeNumber_.at(row);

    uint32_t hashValue;
    switch (hash_) {
        case FilterHash::UNSIGNED_LONG_CAST: {
            hashValue = UnsignedLongIntCast(zm);
            break;
        }

        case FilterHash::BOOST_HASH_COMBINE: {
            hashValue = BoostHashCombine(zm);
            break;
        }

        default:
            throw std::runtime_error{"[pbbam] PBI filter ERROR: unsupported filter hash type"};
    }

    const auto modResult = hashValue % denominator_;
    return Compare::Check(modResult, value_, cmp_);
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PBIFILTERTYPES_INL
