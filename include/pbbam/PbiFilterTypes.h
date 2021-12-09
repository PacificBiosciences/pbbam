#ifndef PBBAM_PBIFILTERTYPES_H
#define PBBAM_PBIFILTERTYPES_H

#include <pbbam/Config.h>

#include <pbbam/Compare.h>
#include <pbbam/PbiFile.h>
#include <pbbam/PbiFilter.h>

#include <boost/optional.hpp>

#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

namespace internal {

/// \internal
///
/// Provides basic container for value/compare-type pair
///
template <typename T>
class FilterBase
{
public:
    T value_;
    boost::optional<std::vector<T>> multiValue_;
    Compare::Type cmp_;

protected:
    FilterBase(T value, const Compare::Type cmp);
    FilterBase(std::vector<T> values, const Compare::Type cmp = Compare::CONTAINS);

    bool CompareHelper(const T& lhs) const;

private:
    bool CompareSingleHelper(const T& lhs) const;
    bool CompareMultiHelper(const T& lhs) const;
};

/// \internal
///
/// Dispatches the lookup to BarcodeLookupData
///
template <typename T, PbiFile::BarcodeField field>
class BarcodeDataFilterBase : public FilterBase<T>
{
public:
    bool Accepts(const PbiRawData& idx, const size_t row) const;

protected:
    BarcodeDataFilterBase(T value, const Compare::Type cmp);
    BarcodeDataFilterBase(std::vector<T> values, const Compare::Type cmp = Compare::CONTAINS);
};

/// \internal
///
/// Dispatches the lookup to BasicLookupData
///
template <typename T, PbiFile::BasicField field>
class BasicDataFilterBase : public FilterBase<T>
{
public:
    bool Accepts(const PbiRawData& idx, const size_t row) const;

protected:
    BasicDataFilterBase(T value, const Compare::Type cmp);
    BasicDataFilterBase(std::vector<T> values, const Compare::Type cmp = Compare::CONTAINS);
};

/// \internal
///
/// Dispatches the lookup to MappedLookupData
///
template <typename T, PbiFile::MappedField field>
class MappedDataFilterBase : public FilterBase<T>
{
public:
    bool Accepts(const PbiRawData& idx, const size_t row) const;

protected:
    MappedDataFilterBase(T value, const Compare::Type cmp);
    MappedDataFilterBase(std::vector<T> values, const Compare::Type cmp = Compare::CONTAINS);
};

}  // namespace internal

/// \brief The PbiAlignedEndFilter class provides a PbiFilter-compatible filter
///        on aligned end.
///
/// Example: \include code/PbiAlignedEndFilter.txt
///
/// \sa BamRecord::AlignedEnd
///
class PbiAlignedEndFilter
    : public internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::A_END>
{
public:
    /// \brief Creates a filter on aligned end.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiAlignedEndFilter(const uint32_t position, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiAlignedLengthFilter class provides a PbiFilter-compatible
///        filter on aligned length.
///
/// Example: \include code/PbiAlignedLengthFilter.txt
///
/// \sa BamRecord::AlignedEnd, BamRecord::AlignedStart
///
class PbiAlignedLengthFilter : public internal::FilterBase<uint32_t>
{
public:
    /// \brief Creates a filter on aligned length.
    ///
    /// \param[in] length value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiAlignedLengthFilter(const uint32_t length, const Compare::Type cmp = Compare::EQUAL);

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
class PbiAlignedStartFilter
    : public internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::A_START>
{
public:
    /// \brief Creates a filter on aligned start.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiAlignedStartFilter(const uint32_t position, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiAlignedStrandFilter class provides a PbiFilter-compatible
///        filter on aligned strand.
///
/// Example: \include code/PbiAlignedStrandFilter.txt
///
/// \sa BamRecord::AlignedStrand
///
class PbiAlignedStrandFilter
    : public internal::MappedDataFilterBase<Data::Strand, PbiFile::MappedField::STRAND>
{
public:
    /// \brief Creates a strand filter.
    ///
    /// \param[in] strand  strand value to compare on
    /// \param[in] cmp     compare type
    ///
    PbiAlignedStrandFilter(const Strand strand, const Compare::Type cmp = Compare::EQUAL);
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
class PbiBarcodeFilter
{
public:
    /// \brief Creates a single-value barcode filter.
    ///
    /// \param[in] barcode  barcode ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeFilter(const int16_t barcode, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted barcode filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no barcodes in the blacklist.
    ///
    /// \param[in] barcodes barcode IDs
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeFilter(std::vector<int16_t> barcodes, const Compare::Type cmp = Compare::CONTAINS);

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
class PbiBarcodeForwardFilter
    : public internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_FORWARD>
{
public:
    /// \brief Creates a single-value forward barcode filter.
    ///
    /// \param[in] bcFwdId  (forward) barcode ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeForwardFilter(const int16_t bcFwdId, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted forward barcode filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] barcodes barcode IDs
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeForwardFilter(std::vector<int16_t> barcodes,
                            const Compare::Type cmp = Compare::CONTAINS);
};

/// \brief The PbiBarcodeQualityFilter class provides a PbiFilter-compatible
///        filter on  barcode quality.
///
/// Example: \include code/PbiBarcodeQualityFilter.txt
///
/// \sa BamRecord::BarcodeQuality
///
class PbiBarcodeQualityFilter
    : public internal::BarcodeDataFilterBase<uint8_t, PbiFile::BarcodeField::BC_QUALITY>
{
public:
    /// \brief Creates a single-value barcode quality filter.
    ///
    /// \param[in] bcQuality    barcode quality to compare on
    /// \param[in] cmp          compare type
    ///
    PbiBarcodeQualityFilter(const uint8_t bcQuality, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiBarcodeReverseFilter class provides a PbiFilter-compatible
///        filter on forward barcode ID.
///
/// Example: \include code/PbiBarcodeReverseFilter.txt
///
/// \sa BamRecord::BarcodeReverse
///
class PbiBarcodeReverseFilter
    : public internal::BarcodeDataFilterBase<int16_t, PbiFile::BarcodeField::BC_REVERSE>
{
public:
    /// \brief Creates a single-value reverse barcode filter.
    ///
    /// \param[in] bcRevId  (reverse) barcode ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeReverseFilter(const int16_t bcRevId, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted reverse barcode filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] barcodes barcode IDs
    /// \param[in] cmp      compare type
    ///
    PbiBarcodeReverseFilter(std::vector<int16_t> barcodes,
                            const Compare::Type cmp = Compare::CONTAINS);
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
class PbiBarcodesFilter
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
    PbiBarcodesFilter(const int16_t bcForward, const int16_t bcReverse,
                      const Compare::Type cmp = Compare::EQUAL);

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
class PbiIdentityFilter : public internal::FilterBase<float>
{
public:
    /// \brief Creates a read identity filter.
    ///
    /// \param[in] identity value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiIdentityFilter(const float identity, const Compare::Type cmp = Compare::EQUAL);

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
class PbiLocalContextFilter
    : public internal::BasicDataFilterBase<Data::LocalContextFlags,
                                           PbiFile::BasicField::CONTEXT_FLAG>
{
public:
    PbiLocalContextFilter(const Data::LocalContextFlags& flags,
                          const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiMapQualityFilter class provides a PbiFilter-compatible filter on
///        mapping quality.
///
/// Example: \include code/PbiMapQualityFilter.txt
///
/// \sa BamRecord::MapQuality
///
class PbiMapQualityFilter
    : public internal::MappedDataFilterBase<uint8_t, PbiFile::MappedField::MAP_QUALITY>
{
public:
    /// \brief Creates a map quality filter.
    ///
    /// \param[in] mapQual  value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiMapQualityFilter(const uint8_t mapQual, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiMovieNameFilter class provides a PbiFilter-compatible filter
///        on movie name.
///
/// Example: \include code/PbiMovieNameFilter.txt
///
/// \sa BamRecord::MovieName
///
class PbiMovieNameFilter
{
public:
    /// \brief Creates a single-value movie name filter.
    ///
    /// \param[in] movieName    movie name to compare on
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match movie name, exactly.
    ///
    PbiMovieNameFilter(const std::string& movieName, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted movie name filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] movieNames   movie names
    /// \param[in] cmp          compare type
    ///
    PbiMovieNameFilter(const std::vector<std::string>& movieNames,
                       const Compare::Type cmp = Compare::CONTAINS);

    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    std::set<std::string> movieNames_;
    mutable std::set<int32_t> candidateRgIds_;
    Compare::Type cmp_;
};

/// \brief The PbiNumDeletedBasesFilter class provides a PbiFilter-compatible
///        filter on the number of deleted bases.
///
/// Example: \include code/PbiNumDeletedBasesFilter.txt
///
/// \sa BamRecord::NumDeletedBases
///
class PbiNumDeletedBasesFilter
    : public internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_DEL>
{
public:
    /// \brief Creates a filter on the number of deleted bases.
    ///
    /// \param[in] numDeletions value to compare on
    /// \param[in] cmp          compare type
    ///
    PbiNumDeletedBasesFilter(const size_t numDeletions, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiNumInsertededBasesFilter class provides a PbiFilter-compatible
///        filter on the number of inserted bases.
///
/// Example: \include code/PbiNumInsertedBasesFilter.txt
///
/// \sa BamRecord::NumInsertedBases
///
class PbiNumInsertedBasesFilter
    : public internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_INS>
{
public:
    /// \brief Creates a filter on the number of inserted bases.
    ///
    /// \param[in] numInsertions    value to compare on
    /// \param[in] cmp              compare type
    ///
    PbiNumInsertedBasesFilter(const size_t numInsertions, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiNumMatchesFilter class provides a PbiFilter-compatible filter
///        on the number of matched bases.
///
/// Example: \include code/PbiNumMatchesFilter.txt
///
/// \sa BamRecord::NumMatches
///
class PbiNumMatchesFilter : public internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_M>
{
public:
    /// \brief Creates a filter on the number of matched bases.
    ///
    /// \param[in] numMatchedBases  value to compare on
    /// \param[in] cmp              compare type
    ///
    PbiNumMatchesFilter(const size_t numMatchedBases, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiNumMismatchesFilter class provides a PbiFilter-compatible
///        filter on the number of mismatched bases.
///
/// Example: \include code/PbiNumMismatchesFilter.txt
///
/// \sa BamRecord::NumMismatches
///
class PbiNumMismatchesFilter
    : public internal::MappedDataFilterBase<size_t, PbiFile::MappedField::N_MM>
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

/// \brief The PbiNumSubreadsFilter class provides a PbiFilter-compatible
///        filter on the number of subreads per ZMW hole number.
///
class PbiNumSubreadsFilter
{
public:
    /// \brief Creates a filter on the number of subreads per ZMW hole number.
    ///
    /// \param[in] numSubreads   value to compare on
    /// \param[in] cmp                  compare type
    ///
    PbiNumSubreadsFilter(int numSubreads, const Compare::Type cmp = Compare::EQUAL);

    PbiNumSubreadsFilter(const PbiNumSubreadsFilter& other);

    PbiNumSubreadsFilter(PbiNumSubreadsFilter&&) noexcept;
    PbiNumSubreadsFilter& operator=(PbiNumSubreadsFilter&&) noexcept;
    ~PbiNumSubreadsFilter();

    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    struct PbiNumSubreadsFilterPrivate;
    std::unique_ptr<PbiNumSubreadsFilterPrivate> d_;
};

/// \brief The PbiQueryEndFilter class provides a PbiFilter-compatible filter
///        on query end.
///
/// Example: \include code/PbiQueryEndFilter.txt
///
/// \sa BamRecord::QueryEnd
///
class PbiQueryEndFilter : public internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::Q_END>
{
public:
    /// \brief Creates a filter on query end position.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiQueryEndFilter(const int32_t position, const Compare::Type cmp = Compare::EQUAL);
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
class PbiQueryLengthFilter : public internal::FilterBase<int32_t>
{
public:
    /// \brief Creates a filter on query length
    ///
    /// \param[in] length   value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiQueryLengthFilter(const int32_t length, const Compare::Type cmp = Compare::EQUAL);

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
class PbiQueryNameFilter
{
public:
    /// \brief Creates a single-value query name filter.
    ///
    /// \param[in] qname    query name to compare on
    ///
    /// \note There is no compare type parameter here, it is always
    ///       Compare::EQUAL. Records will match query name, exactly.
    ///
    PbiQueryNameFilter(const std::string& qname, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted query name filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] queryNames   query names
    /// \param[in] cmp          compare type
    ///
    PbiQueryNameFilter(const std::vector<std::string>& queryNames,
                       const Compare::Type cmp = Compare::CONTAINS);

    PbiQueryNameFilter(const PbiQueryNameFilter& other);

    PbiQueryNameFilter(PbiQueryNameFilter&&) noexcept;
    PbiQueryNameFilter& operator=(PbiQueryNameFilter&&) noexcept;
    ~PbiQueryNameFilter();

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
class PbiQueryStartFilter
    : public internal::BasicDataFilterBase<int32_t, PbiFile::BasicField::Q_START>
{
public:
    /// \brief Creates a filter on query start position.
    ///
    /// \param[in] position value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiQueryStartFilter(const int32_t position, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiReadAccuracyFilter class provides a PbiFilter-compatible filter
///        on read accuracy.
///
/// Example: \include code/PbiReadAccuracyFilter.txt
///
/// \sa BamRecord::ReadAccuracy
///
class PbiReadAccuracyFilter
    : public internal::BasicDataFilterBase<Data::Accuracy, PbiFile::BasicField::READ_QUALITY>
{
public:
    /// \brief Creates a filter on read accuracy.
    ///
    /// \param[in] accuracy value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReadAccuracyFilter(const Data::Accuracy accuracy, const Compare::Type cmp = Compare::EQUAL);
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
class PbiReadGroupFilter
{
public:
    /// \brief Creates a filter on read group (numeric) ID value
    ///
    /// \param[in] rgId     numeric read group ID
    /// \param[in] cmp      compare type
    ///
    /// \sa BamRecord::ReadGroupNumericId
    ///
    PbiReadGroupFilter(const int32_t rgId, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a filter on printable read group ID value
    ///
    /// \param[in] rgId     read group ID string
    /// \param[in] cmp      compare type
    ///
    /// \sa BamRecord::ReadGroupId
    ///
    PbiReadGroupFilter(const std::string& rgId, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a filter on read group (object).
    ///
    /// \param[in] rg   read group object
    /// \param[in] cmp  compare type
    ///
    /// \sa BamRecord::ReadGroup
    ///
    PbiReadGroupFilter(const ReadGroupInfo& rg, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted filter on read group numeric IDs.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] rgIds    numeric read group IDs
    /// \param[in] cmp      compare type
    ///
    PbiReadGroupFilter(const std::vector<int32_t>& rgIds,
                       const Compare::Type cmp = Compare::CONTAINS);

    /// \brief Creates a whitelisted or blacklisted filter on read group string IDs.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] rgIds    read group ID strings
    /// \param[in] cmp      compare type
    ///
    PbiReadGroupFilter(const std::vector<std::string>& rgIds,
                       const Compare::Type cmp = Compare::CONTAINS);

    /// \brief Creates a whitelisted or blacklisted filter on read group objects.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] readGroups   ReadGroupInfo objects
    /// \param[in] cmp          compare type
    ///
    PbiReadGroupFilter(const std::vector<ReadGroupInfo>& readGroups,
                       const Compare::Type cmp = Compare::CONTAINS);

    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    void AddReadGroups(const std::vector<ReadGroupInfo>& readGroups);

    std::unordered_map<int32_t, std::vector<ReadGroupInfo>> readGroups_;
    Compare::Type cmp_;
};

/// \brief The PbiReferenceEndFilter class provides a PbiFilter-compatible
///        filter on reference end.
///
/// Example: \include code/PbiReferenceEndFilter.txt
///
/// \sa BamRecord::ReferenceEnd
///
class PbiReferenceEndFilter
    : public internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::T_END>
{
public:
    /// \brief Creates a filter on reference end.
    ///
    /// \param[in] tEnd     value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReferenceEndFilter(const uint32_t tEnd, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiReferenceIdFilter class provides a PbiFilter-compatible
///        filter on reference ID.
///
/// Example: \include code/PbiReferenceIdFilter.txt
///
/// \sa BamRecord::ReferenceId
///
class PbiReferenceIdFilter
    : public internal::MappedDataFilterBase<int32_t, PbiFile::MappedField::T_ID>
{
public:
    /// \brief Creates a single-value reference ID filter.
    ///
    /// \param[in] tId  reference ID to compare on
    /// \param[in] cmp  compare type
    ///
    PbiReferenceIdFilter(const int32_t tId, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted reference ID filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] tIds reference IDs
    /// \param[in] cmp  compare type
    ///
    PbiReferenceIdFilter(std::vector<int32_t> tIds, const Compare::Type cmp = Compare::CONTAINS);
};

/// \brief The PbiReferenceNameFilter class provides a PbiFilter-compatible
///        filter on reference name.
///
/// Example: \include code/PbiReferenceNameFilter.txt
///
/// \sa BamRecord::ReferenceName
///
class PbiReferenceNameFilter
{
public:
    /// \brief Creates a single-value reference name filter.
    ///
    /// \param[in] rname    reference ID to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReferenceNameFilter(std::string rname, Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted reference name filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] rnames   reference names
    /// \param[in] cmp      compare type
    ///
    PbiReferenceNameFilter(std::vector<std::string> rnames,
                           const Compare::Type cmp = Compare::CONTAINS);

    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    mutable bool initialized_ = false;
    mutable PbiFilter subFilter_;
    std::string rname_;
    boost::optional<std::vector<std::string>> rnameWhitelist_;
    Compare::Type cmp_;

    // marked const so we can delay setup of filter in Accepts(), once we have
    // access to PBI/BAM input. modified values marked mutable accordingly
    void Initialize(const PbiRawData& idx) const;

    void Validate() const;
};

/// \brief The PbiReferenceStartFilter class provides a PbiFilter-compatible
///        filter on reference start.
///
/// Example: \include code/PbiReferenceStartFilter.txt
///
/// \sa BamRecord::ReferenceStart
///
class PbiReferenceStartFilter
    : public internal::MappedDataFilterBase<uint32_t, PbiFile::MappedField::T_START>
{
public:
    /// \brief Creates a filter on reference start.
    ///
    /// \param[in] tStart   value to compare on
    /// \param[in] cmp      compare type
    ///
    PbiReferenceStartFilter(const uint32_t tStart, const Compare::Type cmp = Compare::EQUAL);
};

/// \brief The PbiZmwFilter class provides a PbiFilter-compatible filter on
///        ZMW hole number.
///
/// Example: \include code/PbiZmwFilter.txt
///
/// \sa BamRecord::HoleNumber
///
class PbiZmwFilter
{
public:
    /// \brief Creates a single-value ZMW hole number filter.
    ///
    /// \param[in] zmw  value to compare on
    /// \param[in] cmp  compare type
    ///
    PbiZmwFilter(const int32_t zmw, const Compare::Type cmp = Compare::EQUAL);

    /// \brief Creates a whitelisted or blacklisted ZMW hole number filter.
    ///
    /// \note If compare type is Compare::CONTAINS, accepted records will match
    ///       at least one value from whitelist, in bc_forward.
    ///       If compare type is Compare::NOT_CONTAINS, accepted records will
    ///       match no values in the blacklist.
    ///
    /// \param[in] zmws ZMW hole numbers
    /// \param[in] cmp  compare type
    ///
    PbiZmwFilter(std::vector<int32_t> zmws, const Compare::Type cmp = Compare::CONTAINS);

    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    Compare::Type cmp_;
    int32_t singleZmw_;
    std::set<int32_t> zmwLookup_;
};

// ----------------------------------------------
// NOTE: modulo filtering only enabled for ZMW.
//
// I need to generalize more if we're going to use
// this on more fields.
// ----------------------------------------------

enum class FilterHash
{
    UNSIGNED_LONG_CAST,
    BOOST_HASH_COMBINE,
};

class PbiZmwModuloFilter
{
public:
    PbiZmwModuloFilter(const uint32_t denominator, const uint32_t value,
                       const FilterHash hashtype = FilterHash::UNSIGNED_LONG_CAST,
                       const Compare::Type = Compare::EQUAL);

    /// \brief Performs the actual index lookup.
    ///
    /// Most client code should not need to use this method directly.
    ///
    bool Accepts(const PbiRawData& idx, const size_t row) const;

private:
    uint32_t denominator_;
    uint32_t value_;
    FilterHash hash_;
    Compare::Type cmp_;
};

}  // namespace BAM
}  // namespace PacBio

#include <pbbam/internal/PbiFilterTypes.inl>

#endif  // PBBAM_PBIFILTERTYPES_H
