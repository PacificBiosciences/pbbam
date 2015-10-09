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

// Author: Derek Barnett

#ifndef COMPARE_H
#define COMPARE_H

#include "pbbam/BamRecord.h"
#include <functional>
#include <string>
#include <utility>

namespace PacBio {
namespace BAM {

/// The (namespace-like) struct defines a variety of functor objects that can be
/// used to filter, sort, etc. BamRecord
///
/// \code{.cpp}
///     std::vector<BamRecord> records;
///
///     // sort on increasing ZMW hole number
///     std::sort(records.begin(), records.end(), Compare::Zmw());
///
///     // sort by mapping quality, highest to lowest
///     std::sort(records.begin(), records.end(), Compare::MapQuality(Compare::GREATER_THAN));
///
/// \endcode
///
struct PBBAM_EXPORT Compare
{

public:

    /// \name Comparison Type
    /// \{

    /// Supported comparison types (==, !=, <, <=, >, >=).
    ///
    enum Type {
        EQUAL = 0
      , NOT_EQUAL
      , LESS_THAN
      , LESS_THAN_EQUAL
      , GREATER_THAN
      , GREATER_THAN_EQUAL
    };

    /// Convert operator string to Compare::Type.
    ///
    /// \example
    /// \code {.cpp}
    ///     Compare::Type type = Compare::TypeFromOperator("!=");
    ///     assert(type == Compare::NOT_EQUAL);
    /// \endcode
    ///
    /// \param[in] opString operator string. Can be C++-style operators ("==", "!=", "<=", etc)
    ///            or alpha equivalents ("eq", "ne", "lte", etc).
    /// \returns comparison type from an operator string
    /// \throws std::runtime_error if cannot convert opString to Compare::Type
    /// \sa Compare::TypeToOperator
    ///
    static Compare::Type TypeFromOperator(const std::string& opString);

    /// Convert Compare::Type to printable enum name.
    ///
    /// \example
    /// \code{.cpp}
    ///     string name = Compare::TypeToName(Compare::LESS_THAN);
    ///     assert(name = "Compare::LESS_THAN");
    /// \endcode
    ///
    /// \param[in] type Compare::Type to convert
    /// \returns the printable name for a Compare::Type enum value.are::Type
    /// \throws std::runtime_error on unknown Compare::Type
    ///
    static std::string TypeToName(const Compare::Type& type);

    /// Convert Compare::Type to printable operator.
    ///
    /// \param[in] type Compare::Type to convert
    /// \param[in] asAlpha (optional) flag to print using alpha equivalents e.g. "lte" rather than "<="
    /// \returns the printable operator string
    /// \throws std::runtime_error on unknown Compare::Type
    ///
    static std::string TypeToOperator(const Compare::Type& type, bool asAlpha = false);

    /// \}

public:

    /// \name Comparison Function Objects
    /// \{

    /// Base class for all BamRecord compare functors.
    ///
    /// Mostly used for method signatures that can accept any comparator.
    ///
    /// Custom comparators may be used by inheriting from this class.
    ///
    struct Base : public std::function<bool(const BamRecord&, const BamRecord&)> { };

private:
    /// \internal
    ///
    /// Exists to provide the typedef we'll use in the actual MemberFunctionBase, since we
    /// need to use it in the template signature. This keeps that a lot easier to read.
    ///
    template<typename ValueType>
    struct MemberFunctionBaseHelper : public Compare::Base
    {
        typedef ValueType (BamRecord::*MemberFnType)(void) const;
    };

public:
    /// Base class for all BamRecord compare functors that take a BamRecord function pointer
    /// and compare on its return type.
    ///
    /// Derived comparators usually need only declare the return value & function pointer in
    /// the template signature. This class implements the basic method-calling machinery.
    ///
    /// Custom comparators will work for any BamRecord member function that does not take any input
    /// parameters.
    ///
    template<typename ValueType,
             typename MemberFunctionBaseHelper<ValueType>::MemberFnType fn,
             typename CompareType = std::less<ValueType> >
    struct MemberFunctionBase : public Compare::MemberFunctionBaseHelper<ValueType>
    {
        bool operator()(const BamRecord& lhs, const BamRecord& rhs) const;
    };

public:

    /// Provides an operator() that compares on BamRecord::AlignedEnd
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::AlignedEnd());
    /// \endcode
    ///
    struct AlignedEnd : public MemberFunctionBase<Position, &BamRecord::AlignedEnd> { };

    /// Provides an operator() that compares on BamRecord::AlignedStart
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::AlignedStart());
    /// \endcode
    ///
    struct AlignedStart : public MemberFunctionBase<Position, &BamRecord::AlignedStart> { };

    /// Provides an operator() that compares on BamRecord::AlignedStrand
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::AlignedStrand());
    /// \endcode
    ///
    struct AlignedStrand : public MemberFunctionBase<Strand, &BamRecord::AlignedStrand> { };

    /// Provides an operator() that compares on BamRecord::BarcodeForward
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::BarcodeForward());
    /// \endcode
    ///
    struct BarcodeForward : public MemberFunctionBase<uint16_t, &BamRecord::BarcodeForward> { };

    /// Provides an operator() that compares on BamRecord::BarcodeQuality
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::BarcodeQuality());
    /// \endcode
    ///
    struct BarcodeQuality : public MemberFunctionBase<uint8_t, &BamRecord::BarcodeQuality> { };

    /// Provides an operator() that compares on BamRecord::BarcodeReverse
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::BarcodeReverse());
    /// \endcode
    ///
    struct BarcodeReverse: public MemberFunctionBase<uint16_t, &BamRecord::BarcodeReverse> { };

    /// Provides an operator() that compares on BamRecord::FullName
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::FullName());
    /// \endcode
    ///
    struct FullName : public MemberFunctionBase<std::string, &BamRecord::FullName> { };

    /// Provides an operator() that compares on BamRecord::LocalContextFlags
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::LocalContextFlag());
    /// \endcode
    ///
    struct LocalContextFlag : public MemberFunctionBase<LocalContextFlags, &BamRecord::LocalContextFlags> { };

    /// Provides an operator() that compares on BamRecord::MapQuality
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::MapQuality());
    /// \endcode
    ///
    struct MapQuality : public MemberFunctionBase<uint8_t, &BamRecord::MapQuality> { };

    /// Provides an operator() that compares on BamRecord::MovieName
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::MovieName());
    /// \endcode
    ///
    struct MovieName : public MemberFunctionBase<std::string, &BamRecord::MovieName> { };

    /// Provides an operator() is essentially a no-op for comparing/sorting.
    ///
    /// If used in a sorting operation, then no change will occur.
    ///
    struct None : public Compare::Base
    {
        bool operator()(const BamRecord&, const BamRecord&) const;
    };

    /// Provides an operator() that compares on BamRecord::NumDeletedBases
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::NumDeletedBases());
    /// \endcode
    ///
    struct NumDeletedBases : public MemberFunctionBase<size_t, &BamRecord::NumDeletedBases> { };

    /// Provides an operator() that compares on BamRecord::NumInsertedBases
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::NumInsertedBases());
    /// \endcode
    ///
    struct NumInsertedBases : public MemberFunctionBase<size_t, &BamRecord::NumInsertedBases> { };

    /// Provides an operator() that compares on BamRecord::NumMatches
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::NumMatches());
    /// \endcode
    ///
    struct NumMatches : public MemberFunctionBase<size_t, &BamRecord::NumMatches> { };

    /// Provides an operator() that compares on BamRecord::NumMismatches
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::NumMismatches());
    /// \endcode
    ///
    struct NumMismatches : public MemberFunctionBase<size_t, &BamRecord::NumMismatches> { };

    /// Provides an operator() that compares on BamRecord::QueryEnd
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::QueryEnd());
    /// \endcode
    ///
    struct QueryEnd : public MemberFunctionBase<Position, &BamRecord::QueryEnd> { };

    /// Provides an operator() that compares on BamRecord::QueryStart
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::QueryStart());
    /// \endcode
    ///
    struct QueryStart : public MemberFunctionBase<Position, &BamRecord::QueryStart> { };

    /// Provides an operator() that compares on BamRecord::ReadAccuracy
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::ReadAccuracy());
    /// \endcode
    ///
    struct ReadAccuracy : public MemberFunctionBase<Accuracy, &BamRecord::ReadAccuracy> { };

    /// Provides an operator() that compares on BamRecord::ReadGroupId.
    ///
    /// \note Even though ReadGroupId contains hex values, it is still a std::string.
    ///       Comparisons will use lexical, not numeric ordering. If numeric ordering
    ///       is desired, use Compare::ReadGroupNumericId instead.
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::ReadGroupId());
    /// \endcode
    ///
    struct ReadGroupId : public MemberFunctionBase<std::string, &BamRecord::ReadGroupId> { };

    /// Provides an operator() that compares on BamRecord::ReadGroupNumericId
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::ReadGroupNumericId());
    /// \endcode
    ///
    struct ReadGroupNumericId : public MemberFunctionBase<int32_t, &BamRecord::ReadGroupNumericId> { };

    /// Provides an operator() that compares on BamRecord::ReferenceEnd
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::ReferenceEnd());
    /// \endcode
    ///
    struct ReferenceEnd : public MemberFunctionBase<Position, &BamRecord::ReferenceEnd> { };

    /// Provides an operator() that compares on BamRecord::ReferenceId
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::ReferenceId());
    /// \endcode
    ///
    struct ReferenceId : public MemberFunctionBase<int32_t, &BamRecord::ReferenceId> { };

    /// Provides an operator() that compares on BamRecord::ReferenceName
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::ReferenceName());
    /// \endcode
    ///
    struct ReferenceName : public MemberFunctionBase<std::string, &BamRecord::ReferenceName> { };

    /// Provides an operator() that compares on BamRecord::ReferenceStart
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::ReferenceStart());
    /// \endcode
    ///
    struct ReferenceStart : public MemberFunctionBase<Position, &BamRecord::ReferenceStart> { };

    /// Provides an operator() that compares on BamRecord::HoleNumber
    ///
    /// Example:
    /// \code{.cpp}
    ///     std::vector<BamRecord> records;
    ///     std::sort(records.begin(), records.end(), Compare::Zmw());
    /// \endcode
    ///
    struct Zmw : public MemberFunctionBase<int32_t, &BamRecord::HoleNumber> { };

    /// \}
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/Compare.inl"

#endif // COMPARE_H
