// File Description
/// \file PbiFilter.cpp
/// \brief Implements the PbiFilter class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/PbiFilter.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "StringUtils.h"
#include "pbbam/PbiFilterTypes.h"

namespace PacBio {
namespace BAM {
namespace internal {

// clang-format off
enum class BuiltIn
{
    AlignedEndFilter
  , AlignedLengthFilter
  , AlignedStartFilter
  , AlignedStrandFilter
  , BarcodeFilter
  , BarcodeForwardFilter
  , BarcodeQualityFilter
  , BarcodeReverseFilter
  , BarcodesFilter
  , IdentityFilter
  , LocalContextFilter
  , MovieNameFilter
  , NumDeletedBasesFilter
  , NumInsertedBasesFilter
  , NumMatchesFilter
  , NumMismatchesFilter
  , QueryEndFilter
  , QueryLengthFilter
  , QueryNameFilter
  , QueryNamesFromFileFilter
  , QueryStartFilter
  , ReadAccuracyFilter
  , ReadGroupFilter
  , ReferenceEndFilter
  , ReferenceIdFilter
  , ReferenceNameFilter
  , ReferenceStartFilter
  , ZmwFilter
};

static const std::unordered_map<std::string, BuiltIn> builtInLookup =
{
    // property name   built-in filter
    { "ae",            BuiltIn::AlignedEndFilter },
    { "aend",          BuiltIn::AlignedEndFilter },
    { "alignedlength", BuiltIn::AlignedLengthFilter },
    { "as",            BuiltIn::AlignedStartFilter },
    { "astart",        BuiltIn::AlignedStartFilter },
    { "readstart",     BuiltIn::AlignedStartFilter },
    { "bc",            BuiltIn::BarcodeFilter },
    { "barcode",       BuiltIn::BarcodeFilter },
    { "bcf",           BuiltIn::BarcodeForwardFilter },
    { "bq",            BuiltIn::BarcodeQualityFilter },
    { "bcq",           BuiltIn::BarcodeQualityFilter },
    { "bcr",           BuiltIn::BarcodeReverseFilter },
    { "accuracy",      BuiltIn::IdentityFilter },
    { "identity",      BuiltIn::IdentityFilter },
    { "cx",            BuiltIn::LocalContextFilter },
    { "movie",         BuiltIn::MovieNameFilter },
    { "qe",            BuiltIn::QueryEndFilter },
    { "qend",          BuiltIn::QueryEndFilter },
    { "length",        BuiltIn::QueryLengthFilter },
    { "querylength",   BuiltIn::QueryLengthFilter },
    { "qname",         BuiltIn::QueryNameFilter },
    { "qname_file",    BuiltIn::QueryNamesFromFileFilter },
    { "qs",            BuiltIn::QueryStartFilter },
    { "qstart",        BuiltIn::QueryStartFilter },
    { "rq",            BuiltIn::ReadAccuracyFilter },
    { "te",            BuiltIn::ReferenceEndFilter },
    { "tend",          BuiltIn::ReferenceEndFilter },
    { "rname",         BuiltIn::ReferenceNameFilter },
    { "ts",            BuiltIn::ReferenceStartFilter },
    { "tstart",        BuiltIn::ReferenceStartFilter },
    { "pos",           BuiltIn::ReferenceStartFilter },
    { "zm",            BuiltIn::ZmwFilter },
    { "zmw",           BuiltIn::ZmwFilter }
};

static const std::unordered_map<std::string, LocalContextFlags> contextFlagNames =
{
    { "NO_LOCAL_CONTEXT", LocalContextFlags::NO_LOCAL_CONTEXT },
    { "ADAPTER_BEFORE",   LocalContextFlags::ADAPTER_BEFORE },
    { "ADAPTER_AFTER",    LocalContextFlags::ADAPTER_AFTER },
    { "BARCODE_BEFORE",   LocalContextFlags::BARCODE_BEFORE },
    { "BARCODE_AFTER",    LocalContextFlags::BARCODE_AFTER },
    { "FORWARD_PASS",     LocalContextFlags::FORWARD_PASS },
    { "REVERSE_PASS",     LocalContextFlags::REVERSE_PASS }
};
// clang-format off

// helper methods (for handling maybe-list strings))
static inline bool isBracketed(const std::string& value)
{
    static const std::string openBrackets = "[({";
    static const std::string closeBrackets = "])}";
    return openBrackets.find(value.at(0)) != std::string::npos &&
           closeBrackets.find(value.at(value.length() - 1)) != std::string::npos;
}

static inline bool isList(const std::string& value) { return value.find(',') != std::string::npos; }

static PbiFilter CreateBarcodeFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) throw std::runtime_error{"empty value for barcode filter property"};

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        std::vector<std::string> barcodes = internal::Split(value, ',');
        if (barcodes.size() != 2) throw std::runtime_error{"only 2 barcode values expected"};
        return PbiBarcodesFilter{boost::numeric_cast<int16_t>(std::stoi(barcodes.at(0))),
                                 boost::numeric_cast<int16_t>(std::stoi(barcodes.at(1))),
                                 compareType};
    } else
        return PbiBarcodeFilter{boost::numeric_cast<int16_t>(stoi(value)), compareType};
}

static PbiFilter CreateBarcodeForwardFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) throw std::runtime_error{"empty value for barcode_forward filter property"};

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        std::vector<std::string> tokens = internal::Split(value, ',');
        std::vector<int16_t> barcodes;
        barcodes.reserve(tokens.size());
        for (const auto& t : tokens)
            barcodes.push_back(boost::numeric_cast<int16_t>(stoi(t)));
        return PbiBarcodeForwardFilter{std::move(barcodes)};
    } else
        return PbiBarcodeForwardFilter{boost::numeric_cast<int16_t>(std::stoi(value)), compareType};
}

static PbiFilter CreateBarcodeReverseFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) throw std::runtime_error{"empty value for barcode_reverse filter property"};

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        std::vector<std::string> tokens = internal::Split(value, ',');
        std::vector<int16_t> barcodes;
        barcodes.reserve(tokens.size());
        for (const auto& t : tokens)
            barcodes.push_back(boost::numeric_cast<int16_t>(std::stoi(t)));
        return PbiBarcodeReverseFilter{std::move(barcodes)};
    } else
        return PbiBarcodeReverseFilter{boost::numeric_cast<int16_t>(stoi(value)), compareType};
}

static PbiFilter CreateLocalContextFilter(const std::string& value, const Compare::Type compareType)
{
    if (value.empty()) throw std::runtime_error{"empty value for local context filter property"};

    LocalContextFlags filterValue = LocalContextFlags::NO_LOCAL_CONTEXT;

    // if raw integer
    if (isdigit(value.at(0))) filterValue = static_cast<LocalContextFlags>(stoi(value));

    // else interpret as flag names
    else {
        std::vector<std::string> tokens = internal::Split(value, '|');
        for (std::string& token : tokens) {
            boost::algorithm::trim(token);  // trim whitespace
            filterValue = (filterValue | contextFlagNames.at(token));
        }
    }

    return PbiFilter{PbiLocalContextFilter{filterValue, compareType}};
}

static PbiFilter CreateQueryNamesFilterFromFile(const std::string& value, const DataSet& dataset)
{
    // resolve file from dataset, value
    const std::string resolvedFilename = dataset.ResolvePath(value);
    std::vector<std::string> whitelist;
    std::string fn;
    std::ifstream in(resolvedFilename);
    while (std::getline(in, fn))
        whitelist.push_back(fn);
    return PbiQueryNameFilter{whitelist};
}

static PbiFilter CreateZmwFilter(std::string value, const Compare::Type compareType)
{

    if (value.empty()) throw std::runtime_error{"empty value for ZMW filter property"};

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        std::vector<std::string> tokens = internal::Split(value, ',');
        std::vector<int32_t> zmws;
        zmws.reserve(tokens.size());
        for (const auto& t : tokens)
            zmws.push_back(boost::numeric_cast<int32_t>(stoi(t)));
        return PbiZmwFilter{std::move(zmws)};
    } else
        return PbiZmwFilter{boost::numeric_cast<int32_t>(stoi(value)), compareType};
}

static PbiFilter FromDataSetProperty(const Property& property, const DataSet& dataset)
{
    try {
        const std::string& value = property.Value();
        const Compare::Type compareType = Compare::TypeFromOperator(property.Operator());
        const BuiltIn builtInCode =
            builtInLookup.at(boost::algorithm::to_lower_copy(property.Name()));

        // clang-format off
        switch (builtInCode) {

            // single-value filters
            case BuiltIn::AlignedEndFilter     : return PbiAlignedEndFilter{ static_cast<uint32_t>(std::stoul(value)), compareType };
            case BuiltIn::AlignedLengthFilter  : return PbiAlignedLengthFilter{ static_cast<uint32_t>(std::stoul(value)), compareType };
            case BuiltIn::AlignedStartFilter   : return PbiAlignedStartFilter{ static_cast<uint32_t>(std::stoul(value)), compareType };
            case BuiltIn::BarcodeQualityFilter : return PbiBarcodeQualityFilter{ static_cast<uint8_t>(std::stoul(value)), compareType };
            case BuiltIn::IdentityFilter       : return PbiIdentityFilter{ std::stof(value), compareType };
            case BuiltIn::MovieNameFilter      : return PbiMovieNameFilter{ value };
            case BuiltIn::QueryEndFilter       : return PbiQueryEndFilter{ std::stoi(value), compareType };
            case BuiltIn::QueryLengthFilter    : return PbiQueryLengthFilter{ std::stoi(value), compareType };
            case BuiltIn::QueryNameFilter      : return PbiQueryNameFilter{ value };
            case BuiltIn::QueryStartFilter     : return PbiQueryStartFilter{ std::stoi(value), compareType };
            case BuiltIn::ReadAccuracyFilter   : return PbiReadAccuracyFilter{ std::stof(value), compareType };
            case BuiltIn::ReadGroupFilter      : return PbiReadGroupFilter{ value, compareType };
            case BuiltIn::ReferenceEndFilter   : return PbiReferenceEndFilter{ static_cast<uint32_t>(std::stoul(value)), compareType };
            case BuiltIn::ReferenceIdFilter    : return PbiReferenceIdFilter{ std::stoi(value), compareType };
            case BuiltIn::ReferenceNameFilter  : return PbiReferenceNameFilter{ value };
            case BuiltIn::ReferenceStartFilter : return PbiReferenceStartFilter{ static_cast<uint32_t>(std::stoul(value)), compareType };

            // (maybe) list-value filters
            case BuiltIn::BarcodeFilter        : return CreateBarcodeFilter(value, compareType);
            case BuiltIn::BarcodeForwardFilter : return CreateBarcodeForwardFilter(value, compareType);
            case BuiltIn::BarcodeReverseFilter : return CreateBarcodeReverseFilter(value, compareType);
            case BuiltIn::LocalContextFilter   : return CreateLocalContextFilter(value, compareType);
            case BuiltIn::ZmwFilter            : return CreateZmwFilter(value, compareType);

            // other built-ins
            case BuiltIn::QueryNamesFromFileFilter : return CreateQueryNamesFilterFromFile(value, dataset); // compareType ignored

            default :
            throw std::runtime_error{""};
        }
        // clang-format on

        // unreachable
        return PbiFilter{};

    } catch (std::exception& e) {
        std::ostringstream s;
        s << "error: could not create filter from XML Property element:\n"
          << "  Name:     " << property.Name() << '\n'
          << "  Value:    " << property.Value() << '\n'
          << "  Operator: " << property.Operator() << '\n'
          << "  reason:   " << e.what() << '\n';
        throw std::runtime_error{s.str()};
    }
}

}  // namespace internal

PbiFilter PbiFilter::FromDataSet(const DataSet& dataset)
{
    PbiFilter datasetFilter{PbiFilter::UNION};
    for (const auto& xmlFilter : dataset.Filters()) {
        PbiFilter propertiesFilter;
        for (const auto& xmlProperty : xmlFilter.Properties())
            propertiesFilter.Add(internal::FromDataSetProperty(xmlProperty, dataset));
        datasetFilter.Add(propertiesFilter);
    }
    return datasetFilter;
}

PbiFilter PbiFilter::Intersection(std::vector<PbiFilter> filters)
{
    auto result = PbiFilter{PbiFilter::INTERSECT};
    result.Add(std::move(filters));
    return result;
}

PbiFilter PbiFilter::Union(std::vector<PbiFilter> filters)
{
    auto result = PbiFilter{PbiFilter::UNION};
    result.Add(std::move(filters));
    return result;
}

}  // namespace BAM
}  // namespace PacBio
