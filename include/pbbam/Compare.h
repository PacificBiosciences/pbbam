#ifndef PBBAM_COMPARE_H
#define PBBAM_COMPARE_H

#include <pbbam/Config.h>

#include <pbbam/BamRecord.h>

#include <functional>
#include <string>
#include <utility>

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief The Compare class provides utilities for sorting collections of
///        BamRecords.
///
/// \note The functors provided here currently only support std::less<T>
///       comparisons (i.e. sorting by ascending value).
///
/// \include code/Compare.txt
///
struct PBBAM_EXPORT Compare
{
public:
    /// \name Comparison Type
    /// \{

    /// \brief This enum defines the supported comparison types
    ///        { ==, !=, <, <=, >, >=, & (contains), ~ (not contains) }.
    ///
    enum Type
    {
        EQUAL = 0,
        NOT_EQUAL,
        LESS_THAN,
        LESS_THAN_EQUAL,
        GREATER_THAN,
        GREATER_THAN_EQUAL,
        CONTAINS,
        NOT_CONTAINS
    };

    /// \brief Convert operator string to Compare::Type.
    ///
    /// \include code/Compare_TypeFromOperator.txt
    ///
    /// \param[in] opString operator string. Can be C++-style operators
    ///                     ("==", "!=", "<=", etc) or alpha equivalents
    ///                     ("eq", "ne", "lte", etc).
    ///
    /// \returns comparison type from an operator string
    /// \throws std::runtime_error if cannot convert opString to Compare::Type
    /// \sa Compare::TypeToOperator
    ///
    static Compare::Type TypeFromOperator(const std::string& opString);

    /// \brief Convert a Compare::Type to printable enum name.
    ///
    /// \include code/Compare_TypeToName.txt
    ///
    /// \param[in] type Compare::Type to convert
    /// \returns the printable name for a Compare::Type enum value.are::Type
    /// \throws std::runtime_error on unknown Compare::Type
    ///
    static std::string TypeToName(const Compare::Type& type);

    /// \brief Convert a Compare::Type to printable operator.
    ///
    /// \param[in] type     Compare::Type to convert
    /// \param[in] asAlpha  (optional) flag to print using alpha equivalents
    ///                     e.g. "lte" rather than "<="
    /// \returns the printable operator string
    /// \throws std::runtime_error on unknown Compare::Type
    ///
    static std::string TypeToOperator(const Compare::Type& type, bool asAlpha = false);

    /// \}

public:
    /// \name Comparison Function Objects
    /// \{

    /// %Base class for all BamRecord compare functors.
    ///
    /// Mostly used for method signatures that can accept any comparator.
    ///
    /// Custom comparators may be used by inheriting from this class.
    ///
    struct Base : public std::function<bool(const BamRecord&, const BamRecord&)>
    {};

private:
    /// \internal
    ///
    /// Exists to provide the typedef we'll use in the actual
    /// MemberFunctionBase, since we need to use it in the template signature.
    /// This keeps that a lot easier to read.
    ///
    template <typename ValueType>
    struct MemberFunctionBaseHelper : public Compare::Base
    {
        using MemberFnType = ValueType (BamRecord::*)() const;
    };

public:
    /// \brief %Base class for all BamRecord compare functors that take a
    ///        BamRecord function pointer and compare on its return type.
    ///
    /// Derived comparators usually need only declare the return value &
    /// function pointer in the template signature. This class implements the
    /// basic method-calling machinery.
    ///
    /// Custom comparators will work for any BamRecord member function that does
    /// not take any input parameters.
    ///
    template <typename ValueType, typename MemberFunctionBaseHelper<ValueType>::MemberFnType fn,
              typename CompareType = std::less<ValueType> >
    struct MemberFunctionBase : public Compare::MemberFunctionBaseHelper<ValueType>
    {
        bool operator()(const BamRecord& lhs, const BamRecord& rhs) const;
    };

public:
    /// \brief Compares on BamRecord::AlignedEnd.
    ///
    /// Example:
    /// \include code/Compare_AlignedEnd.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct AlignedEnd : public MemberFunctionBase<Data::Position, &BamRecord::AlignedEnd>
    {};

    /// \brief Compares on BamRecord::AlignedStart.
    ///
    /// Example:
    /// \include code/Compare_AlignedStart.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct AlignedStart : public MemberFunctionBase<Data::Position, &BamRecord::AlignedStart>
    {};

    /// \brief Compares on BamRecord::AlignedStrand
    ///
    /// Example:
    /// \include code/Compare_AlignedStrand.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct AlignedStrand : public MemberFunctionBase<Data::Strand, &BamRecord::AlignedStrand>
    {};

    /// \brief Compares on reference ID, then by position.
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct AlignmentPosition : public Compare::Base
    {
        bool operator()(const BamRecord& lhs, const BamRecord& rhs) const;
    };

    /// \brief Compares on BamRecord::BarcodeForward.
    ///
    /// Example:
    /// \include code/Compare_BarcodeForward.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct BarcodeForward : public MemberFunctionBase<std::int16_t, &BamRecord::BarcodeForward>
    {};

    /// \brief Compares on BamRecord::BarcodeQuality.
    ///
    /// Example:
    /// \include code/Compare_BarcodeQuality.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct BarcodeQuality : public MemberFunctionBase<std::uint8_t, &BamRecord::BarcodeQuality>
    {};

    /// \brief Compares on BamRecord::BarcodeReverse.
    ///
    /// Example:
    /// \include code/Compare_BarcodeReverse.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct BarcodeReverse : public MemberFunctionBase<std::int16_t, &BamRecord::BarcodeReverse>
    {};

    /// \brief Compares on BamRecord::FullName (lexicographical).
    ///
    /// For PacBio BAM standard-aware sorting on QNAME, use Compare::QName.
    ///
    /// Example:
    /// \include code/Compare_FullName.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct FullName : public MemberFunctionBase<std::string, &BamRecord::FullName>
    {};

    /// \brief Compares on BamRecord::LocalContextFlags.
    ///
    /// Example:
    /// \include code/Compare_LocalContextFlag.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct LocalContextFlag
        : public MemberFunctionBase<Data::LocalContextFlags, &BamRecord::LocalContextFlags>
    {};

    /// \brief Compares on BamRecord::MapQuality.
    ///
    /// Example:
    /// \include code/Compare_MapQuality.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct MapQuality : public MemberFunctionBase<std::uint8_t, &BamRecord::MapQuality>
    {};

    /// \brief Compares on BamRecord::MovieName.
    ///
    /// Example:
    /// \include code/Compare_MovieName.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct MovieName : public MemberFunctionBase<std::string, &BamRecord::MovieName>
    {};

    /// \brief Provides an operator() is essentially a no-op for
    ///        comparing/sorting.
    ///
    /// If used in a sorting operation, then no change will occur.
    ///
    struct None : public Compare::Base
    {
        bool operator()(const BamRecord&, const BamRecord&) const noexcept;
    };

    ///\brief Compares on BamRecord::NumDeletedBases.
    ///
    /// Example:
    /// \include code/Compare_NumDeletedBases.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct NumDeletedBases : public MemberFunctionBase<std::size_t, &BamRecord::NumDeletedBases>
    {};

    /// \brief Compares on BamRecord::NumInsertedBases.
    ///
    /// Example:
    /// \include code/Compare_NumInsertedBases.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct NumInsertedBases : public MemberFunctionBase<std::size_t, &BamRecord::NumInsertedBases>
    {};

    /// \brief Compares on BamRecord::NumMatches.
    ///
    /// Example:
    /// \include code/Compare_NumMatches.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct NumMatches : public MemberFunctionBase<std::size_t, &BamRecord::NumMatches>
    {};

    /// \brief Compares on BamRecord::NumMismatches.
    ///
    /// Example:
    /// \include code/Compare_NumMismatches.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct NumMismatches : public MemberFunctionBase<std::size_t, &BamRecord::NumMismatches>
    {};

    /// \brief Compares BamRecords' QNAMEs, via PacBio BAM spec-aware
    ///        sorting order.
    ///
    /// For lexicographical sorting on QNAME, use Compare::FullName.
    ///
    /// \note Only supports sorting by ascending value, per the PacBio BAM spec.
    ///
    struct QName : public Compare::Base
    {
        bool operator()(const BamRecord& lhs, const BamRecord& rhs) const;
    };

    /// \brief Compares on BamRecord::QueryEnd.
    ///
    /// Example:
    /// \include code/Compare_QueryEnd.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct QueryEnd : public MemberFunctionBase<Data::Position, &BamRecord::QueryEnd>
    {};

    /// \brief Compares on BamRecord::QueryStart.
    ///
    /// Example:
    /// \include code/Compare_QueryStart.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct QueryStart : public MemberFunctionBase<Data::Position, &BamRecord::QueryStart>
    {};

    /// \brief Compares on BamRecord::ReadAccuracy.
    ///
    /// Example:
    /// \include code/Compare_ReadAccuracy.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct ReadAccuracy : public MemberFunctionBase<Data::Accuracy, &BamRecord::ReadAccuracy>
    {};

    /// \brief Compares on BamRecord::ReadGroupId.
    ///
    /// \note Even though the ReadGroupId string contains hex values, it is
    ///       still just a std::string. Comparisons will use lexical, not
    ///       numeric ordering. If numeric ordering is desired, use
    ///       Compare::ReadGroupNumericId instead.
    ///
    /// Example:
    /// \include code/Compare_ReadGroupId.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct ReadGroupId : public MemberFunctionBase<std::string, &BamRecord::ReadGroupId>
    {};

    /// \brief Compares on BamRecord::ReadGroupNumericId.
    ///
    /// Example:
    /// \include code/Compare_ReadGroupNumericId.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct ReadGroupNumericId
        : public MemberFunctionBase<std::int32_t, &BamRecord::ReadGroupNumericId>
    {};

    /// \brief Compares on BamRecord::ReferenceEnd.
    ///
    /// Example:
    /// \include code/Compare_ReferenceEnd.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct ReferenceEnd : public MemberFunctionBase<Data::Position, &BamRecord::ReferenceEnd>
    {};

    /// \brief Compares on BamRecord::ReferenceId.
    ///
    /// Example:
    /// \include code/Compare_ReferenceId.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct ReferenceId : public MemberFunctionBase<std::int32_t, &BamRecord::ReferenceId>
    {};

    /// \brief Compares on BamRecord::ReferenceName.
    ///
    /// Example:
    /// \include code/Compare_ReferenceName.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct ReferenceName : public MemberFunctionBase<std::string, &BamRecord::ReferenceName>
    {};

    /// \brief Compares on BamRecord::ReferenceStart.
    ///
    /// Example:
    /// \include code/Compare_ReferenceStart.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct ReferenceStart : public MemberFunctionBase<Data::Position, &BamRecord::ReferenceStart>
    {};

    /// \brief Compares on BamRecord::HoleNumber.
    ///
    /// Example:
    /// \include code/Compare_Zmw.txt
    ///
    /// \note Currently only supports std::less<T> comparisons (i.e. sorting by
    ///       ascending value).
    ///
    struct Zmw : public MemberFunctionBase<std::int32_t, &BamRecord::HoleNumber>
    {};

    /// \}

    template <typename T>
    static bool Check(const T& lhs, const T& rhs, const Compare::Type cmp)
    {
        switch (cmp) {
            case Compare::EQUAL:
            case Compare::CONTAINS:
                return lhs == rhs;
            case Compare::LESS_THAN:
                return lhs < rhs;
            case Compare::LESS_THAN_EQUAL:
                return lhs <= rhs;
            case Compare::GREATER_THAN:
                return lhs > rhs;
            case Compare::GREATER_THAN_EQUAL:
                return lhs >= rhs;
            case Compare::NOT_EQUAL:
            case Compare::NOT_CONTAINS:
                return lhs != rhs;
            default:
                assert(false);
                throw std::runtime_error{
                    "[pbbam] compare ERROR: encountered unsupported compare type"};
        }
    }
};

}  // namespace BAM
}  // namespace PacBio

#include <pbbam/internal/Compare.inl>

#endif  // PBBAM_COMPARE_H
