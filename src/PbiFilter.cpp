#include "PbbamInternalConfig.h"

#include <pbbam/PbiFilter.h>

#include <pbbam/PbiFilterTypes.h>
#include <pbbam/StringUtilities.h>
#include "FileUtils.h"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include <cctype>
#include <cstdint>

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
  , MapQualityFilter
  , MovieNameFilter
  , NumDeletedBasesFilter
  , NumInsertedBasesFilter
  , NumMatchesFilter
  , NumMismatchesFilter
  , NumSubreadsFilter
  , QIdFilter
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
    { "mapqv",         BuiltIn::MapQualityFilter },
    { "movie",         BuiltIn::MovieNameFilter },
    { "n_subreads",    BuiltIn::NumSubreadsFilter },
    { "qid",           BuiltIn::QIdFilter },
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

static const std::unordered_map<std::string, Data::LocalContextFlags> contextFlagNames =
{
    { "NO_LOCAL_CONTEXT",   Data::LocalContextFlags::NO_LOCAL_CONTEXT },
    { "ADAPTER_BEFORE",     Data::LocalContextFlags::ADAPTER_BEFORE },
    { "ADAPTER_AFTER",      Data::LocalContextFlags::ADAPTER_AFTER },
    { "BARCODE_BEFORE",     Data::LocalContextFlags::BARCODE_BEFORE },
    { "BARCODE_AFTER",      Data::LocalContextFlags::BARCODE_AFTER },
    { "FORWARD_PASS",       Data::LocalContextFlags::FORWARD_PASS },
    { "REVERSE_PASS",       Data::LocalContextFlags::REVERSE_PASS },
    { "ADAPTER_BEFORE_BAD", Data::LocalContextFlags::ADAPTER_BEFORE_BAD},
    { "ADAPTER_AFTER_BAD",  Data::LocalContextFlags::ADAPTER_AFTER_BAD}
};
// clang-format on

// helper methods (for handling maybe-list strings))
static bool isBracketed(const std::string& value)
{
    static const std::string openBrackets{"[({"};
    static const std::string closeBrackets{"])}"};
    return openBrackets.find(value.at(0)) != std::string::npos &&
           closeBrackets.find(value.at(value.length() - 1)) != std::string::npos;
}

static bool isList(const std::string& value) { return value.find(',') != std::string::npos; }

static PbiFilter CreateBarcodeFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: empty value for barcode filter property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        std::vector<std::string> barcodes = Split(value, ',');
        if (barcodes.size() != 2) {
            throw std::runtime_error{"[pbbam] PBI filter ERROR: only 2 barcode values expected"};
        }
        return PbiBarcodesFilter{boost::numeric_cast<int16_t>(std::stoi(barcodes.at(0))),
                                 boost::numeric_cast<int16_t>(std::stoi(barcodes.at(1))),
                                 compareType};
    } else {
        return PbiBarcodeFilter{boost::numeric_cast<int16_t>(std::stoi(value)), compareType};
    }
}

static PbiFilter CreateBarcodeForwardFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: empty value for barcode_forward filter property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        std::vector<std::string> tokens = Split(value, ',');
        std::vector<int16_t> barcodes;
        barcodes.reserve(tokens.size());
        for (const auto& t : tokens) {
            barcodes.push_back(boost::numeric_cast<int16_t>(std::stoi(t)));
        }
        return PbiBarcodeForwardFilter{std::move(barcodes)};
    } else {
        return PbiBarcodeForwardFilter{boost::numeric_cast<int16_t>(std::stoi(value)), compareType};
    }
}

static PbiFilter CreateBarcodeReverseFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: empty value for barcode_reverse filter property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        std::vector<std::string> tokens = Split(value, ',');
        std::vector<int16_t> barcodes;
        barcodes.reserve(tokens.size());
        for (const auto& t : tokens) {
            barcodes.push_back(boost::numeric_cast<int16_t>(std::stoi(t)));
        }
        return PbiBarcodeReverseFilter{std::move(barcodes)};
    } else {
        return PbiBarcodeReverseFilter{boost::numeric_cast<int16_t>(std::stoi(value)), compareType};
    }
}

static PbiFilter CreateLocalContextFilter(const std::string& value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: empty value for local context filter property"};
    }

    Data::LocalContextFlags filterValue = Data::LocalContextFlags::NO_LOCAL_CONTEXT;

    // if raw integer
    if (isdigit(value.at(0))) {
        filterValue = static_cast<Data::LocalContextFlags>(std::stoi(value));

        // else interpret as flag names
    } else {
        std::vector<std::string> tokens = Split(value, '|');
        for (std::string& token : tokens) {
            boost::algorithm::trim(token);  // trim whitespace
            filterValue = (filterValue | contextFlagNames.at(token));
        }
    }

    return PbiFilter{PbiLocalContextFilter{filterValue, compareType}};
}

static PbiFilter CreateMovieNameFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: empty value for movie property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        if (compareType != Compare::EQUAL && compareType != Compare::NOT_EQUAL) {
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unsupported compare type on movie property"};
        }

        std::vector<std::string> tokens = Split(value, ',');
        return PbiMovieNameFilter{std::move(tokens), compareType};
    } else {
        return PbiMovieNameFilter{value, compareType};
    }
}

static PbiFilter CreateQIdFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: empty value for qid property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        if (compareType != Compare::EQUAL && compareType != Compare::NOT_EQUAL) {
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unsupported compare type on qid property"};
        }

        std::vector<int32_t> rgIds;
        for (const auto& t : Split(value, ',')) {
            rgIds.push_back(std::stoi(t));
        }
        return PbiReadGroupFilter{rgIds, compareType};
    } else {
        const int32_t n = std::stoi(value);
        return PbiReadGroupFilter{n, compareType};
    }
}

static PbiFilter CreateQueryNamesFilterFromFile(const std::string& value, const DataSet& dataset,
                                                const Compare::Type compareType)
{
    if (compareType != Compare::EQUAL && compareType != Compare::NOT_EQUAL) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: unsupported compare type on query name property"};
    }

    // resolve file from dataset, value
    const auto resolvedFilename = dataset.ResolvePath(value);
    std::vector<std::string> whitelist;
    std::string fn;
    std::ifstream in(resolvedFilename);
    while (std::getline(in, fn)) {
        whitelist.push_back(fn);
    }
    return PbiQueryNameFilter{whitelist, compareType};
}

static PbiFilter CreateQueryNameFilter(std::string value, const DataSet& dataset,
                                       const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: empty value for query name property"};
    }

    // try possible filename first
    const auto resolvedFilename = dataset.ResolvePath(value);
    if (FileUtils::Exists(value)) {
        return CreateQueryNamesFilterFromFile(value, dataset, compareType);
    }

    // otherwise "normal" qname (single, or list)
    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        if (compareType != Compare::EQUAL && compareType != Compare::NOT_EQUAL) {
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unsupported compare type on query name property"};
        }

        std::vector<std::string> tokens = Split(value, ',');
        return PbiQueryNameFilter{std::move(tokens), compareType};
    } else {
        return PbiQueryNameFilter{value, compareType};
    }
}

static PbiFilter CreateReadGroupFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: empty value for read group property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        if (compareType != Compare::EQUAL && compareType != Compare::NOT_EQUAL) {
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unsupported compare type on read group property"};
        }

        std::vector<std::string> tokens = Split(value, ',');
        return PbiReadGroupFilter{std::move(tokens), compareType};
    } else {
        return PbiReadGroupFilter{value, compareType};
    }
}

static PbiFilter CreateReferenceIdFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: empty value for reference ID property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        if (compareType != Compare::EQUAL && compareType != Compare::NOT_EQUAL) {
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unsupported compare type on reference name ID property"};
        }

        const std::vector<std::string> tokens = Split(value, ',');
        std::vector<int32_t> ids;
        ids.reserve(tokens.size());
        for (const auto& t : tokens) {
            ids.push_back(boost::numeric_cast<int32_t>(std::stoi(t)));
        }
        return PbiReferenceIdFilter{std::move(ids), compareType};
    } else {
        return PbiReferenceIdFilter{boost::numeric_cast<int32_t>(std::stoi(value)), compareType};
    }
}

static PbiFilter CreateReferenceNameFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: empty value for reference name property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        if (compareType != Compare::EQUAL && compareType != Compare::NOT_EQUAL) {
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unsupported compare type on reference name property"};
        }

        std::vector<std::string> tokens = Split(value, ',');
        return PbiReferenceNameFilter{std::move(tokens), compareType};
    } else {
        return PbiReferenceNameFilter{value, compareType};
    }
}

static PbiFilter CreateZmwFilter(std::string value, const Compare::Type compareType)
{
    if (value.empty()) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: empty value for ZMW filter property"};
    }

    if (isBracketed(value)) {
        value.erase(0, 1);
        value.pop_back();
    }

    if (isList(value)) {
        const std::vector<std::string> tokens = Split(value, ',');
        std::vector<int32_t> zmws;
        zmws.reserve(tokens.size());
        for (const auto& t : tokens) {
            zmws.push_back(boost::numeric_cast<int32_t>(std::stoi(t)));
        }
        return PbiZmwFilter{std::move(zmws)};
    } else {
        return PbiZmwFilter{boost::numeric_cast<int32_t>(std::stoi(value)), compareType};
    }
}

static PbiFilter CreateZmwModuloFilter(const Property& property)
{
    if (!property.HasAttribute("Modulo") || !property.HasAttribute("Hash") ||
        property.Name() != "zm") {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: modulo filter is not supported on property: " +
            property.Name()};
    }

    const auto hashType = property.Attribute("Hash");
    const FilterHash hash = [&hashType]() {
        if (boost::algorithm::to_lower_copy(hashType) == "uint32cast") {
            return FilterHash::UNSIGNED_LONG_CAST;
        }
        if (boost::algorithm::to_lower_copy(hashType) == "boosthashcombine") {
            return FilterHash::BOOST_HASH_COMBINE;
        }
        throw std::runtime_error{"[pbbam] PBI filter ERROR: unsupported hash type: " + hashType};
    }();

    const uint32_t denom = std::stoul(property.Attribute("Modulo"));
    const uint32_t value = std::stoul(property.Value());

    return PbiZmwModuloFilter{denom, value, hash, Compare::EQUAL};
}

static PbiFilter FromDataSetProperty(const Property& property, const DataSet& dataset)
{
    try {
        const std::string& value = property.Value();
        if (property.Name() == "zm" && property.HasAttribute("Modulo")) {
            return CreateZmwModuloFilter(property);
        }

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
            case BuiltIn::MapQualityFilter     : return PbiMapQualityFilter{ static_cast<uint8_t>(std::stoul(value)), compareType };
            case BuiltIn::NumSubreadsFilter    : return PbiNumSubreadsFilter{ std::stoi(value), compareType };
            case BuiltIn::QueryEndFilter       : return PbiQueryEndFilter{ std::stoi(value), compareType };
            case BuiltIn::QueryLengthFilter    : return PbiQueryLengthFilter{ std::stoi(value), compareType };
            case BuiltIn::QueryStartFilter     : return PbiQueryStartFilter{ std::stoi(value), compareType };
            case BuiltIn::ReadAccuracyFilter   : return PbiReadAccuracyFilter{ std::stof(value), compareType };
            case BuiltIn::ReferenceEndFilter   : return PbiReferenceEndFilter{ static_cast<uint32_t>(std::stoul(value)), compareType };
            case BuiltIn::ReferenceStartFilter : return PbiReferenceStartFilter{ static_cast<uint32_t>(std::stoul(value)), compareType };

            // (maybe) list-value filters
            case BuiltIn::BarcodeFilter        : return CreateBarcodeFilter(value, compareType);
            case BuiltIn::BarcodeForwardFilter : return CreateBarcodeForwardFilter(value, compareType);
            case BuiltIn::BarcodeReverseFilter : return CreateBarcodeReverseFilter(value, compareType);
            case BuiltIn::LocalContextFilter   : return CreateLocalContextFilter(value, compareType);
            case BuiltIn::MovieNameFilter      : return CreateMovieNameFilter(value, compareType);
            case BuiltIn::QIdFilter            : return CreateQIdFilter(value, compareType);
            case BuiltIn::QueryNameFilter      : return CreateQueryNameFilter(value, dataset, compareType);
            case BuiltIn::ReadGroupFilter      : return CreateReadGroupFilter(value, compareType);
            case BuiltIn::ReferenceIdFilter    : return CreateReferenceIdFilter(value, compareType);
            case BuiltIn::ReferenceNameFilter  : return CreateReferenceNameFilter(value, compareType);
            case BuiltIn::ZmwFilter            : return CreateZmwFilter(value, compareType);

            // other built-ins
            case BuiltIn::QueryNamesFromFileFilter : return CreateQueryNamesFilterFromFile(value, dataset, compareType);

            default :
            throw std::runtime_error{"[pbbam] PBI filter ERROR: invalid built-in filter requested"};
        }
        // clang-format on

        // unreachable
        return PbiFilter{};

    } catch (std::exception& e) {
        std::ostringstream s;
        s << "[pbbam] PBI filter ERROR: could not create filter from XML Property element:\n"
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
        for (const auto& xmlProperty : xmlFilter.Properties()) {
            propertiesFilter.Add(internal::FromDataSetProperty(xmlProperty, dataset));
        }
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
