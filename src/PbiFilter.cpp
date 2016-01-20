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
/// \file PbiFilter.cpp
/// \brief Implements the PbiFilter class.
//
// Author: Derek Barnett

#include "pbbam/PbiFilter.h"
#include "pbbam/PbiFilterTypes.h"
#include <boost/algorithm/string/case_conv.hpp>
#include <algorithm>
#include <sstream>
#include <string>
#include <unordered_map>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

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
  , MovieNameFilter
  , NumDeletedBasesFilter
  , NumInsertedBasesFilter
  , NumMatchesFilter
  , NumMismatchesFilter
  , QueryEndFilter
  , QueryLengthFilter
  , QueryNameFilter
  , QueryStartFilter
  , ReadAccuracyFilter
  , ReadGroupFilter
  , ReferenceEndFilter
  , ReferenceIdFilter
  , ReferenceNameFilter
  , ReferenceStartFilter
  , ZmwFilter
};

static const unordered_map<string, BuiltIn> builtInLookup =
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
    { "accuracy",      BuiltIn::IdentityFilter },
    { "identity",      BuiltIn::IdentityFilter },
    { "movie",         BuiltIn::MovieNameFilter },
    { "qe",            BuiltIn::QueryEndFilter },
    { "qend",          BuiltIn::QueryEndFilter },
    { "length",        BuiltIn::QueryLengthFilter },
    { "querylength",   BuiltIn::QueryLengthFilter },
    { "qname",         BuiltIn::QueryNameFilter },
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

    // no implementation yet for these below
};

static
PbiFilter FromDataSetProperty(const Property& property)
{
    try {
        const string& value = property.Value();
        const Compare::Type compareType = Compare::TypeFromOperator(property.Operator());
        const BuiltIn builtInCode = builtInLookup.at(boost::algorithm::to_lower_copy(property.Name()));
        switch (builtInCode) {
            case BuiltIn::AlignedEndFilter       : return PbiAlignedEndFilter{ static_cast<uint32_t>(stoul(value)), compareType };
            case BuiltIn::AlignedLengthFilter    : return PbiAlignedLengthFilter{ static_cast<uint32_t>(stoul(value)), compareType };
            case BuiltIn::AlignedStartFilter     : return PbiAlignedStartFilter{ static_cast<uint32_t>(stoul(value)), compareType };
            case BuiltIn::BarcodeFilter          : return PbiBarcodeFilter{ static_cast<uint16_t>(stoul(value)), compareType };
            case BuiltIn::IdentityFilter         : return PbiIdentityFilter{ stof(value), compareType };
            case BuiltIn::MovieNameFilter        : return PbiMovieNameFilter{ value };
            case BuiltIn::QueryEndFilter         : return PbiQueryEndFilter{ stoi(value), compareType };
            case BuiltIn::QueryLengthFilter      : return PbiQueryLengthFilter{ stoi(value), compareType };
            case BuiltIn::QueryNameFilter        : return PbiQueryNameFilter{ value };
            case BuiltIn::QueryStartFilter       : return PbiQueryStartFilter{ stoi(value), compareType };
            case BuiltIn::ReadAccuracyFilter     : return PbiReadAccuracyFilter{ stof(value), compareType };
            case BuiltIn::ReadGroupFilter        : return PbiReadGroupFilter{ value, compareType };
            case BuiltIn::ReferenceEndFilter     : return PbiReferenceEndFilter{ static_cast<uint32_t>(stoul(value)), compareType };
            case BuiltIn::ReferenceIdFilter      : return PbiReferenceIdFilter{ stoi(value), compareType };
            case BuiltIn::ReferenceNameFilter    : return PbiReferenceNameFilter{ value };
            case BuiltIn::ReferenceStartFilter   : return PbiReferenceStartFilter { static_cast<uint32_t>(stoul(value)), compareType };
            case BuiltIn::ZmwFilter              : return PbiZmwFilter { stoi(value), compareType };
            default :
                throw std::exception();
        }
        // unreachable
        return PbiFilter{ };

    } catch (std::exception& e) {
        stringstream s;
        s << "error: could not create filter from XML Property element: " << endl
          << "  Name:     " << property.Name()     << endl
          << "  Value:    " << property.Value()    << endl
          << "  Operator: " << property.Operator() << endl
          << "  reason:   " << e.what() << endl;
        throw std::runtime_error(s.str());
    }
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

PbiFilter PbiFilter::FromDataSet(const DataSet& dataset)
{
    auto datasetFilter = PbiFilter{ PbiFilter::UNION };
    for (auto&& xmlFilter : dataset.Filters()) {
        auto propertiesFilter = PbiFilter{ };
        for (auto&& xmlProperty : xmlFilter.Properties())
            propertiesFilter.Add(internal::FromDataSetProperty(xmlProperty));
        datasetFilter.Add(propertiesFilter);
    }
    return datasetFilter;
}

PbiFilter PbiFilter::Intersection(const std::vector<PbiFilter>& filters)
{
    auto result = PbiFilter{ PbiFilter::INTERSECT };
    result.Add(filters);
    return result;
}

PbiFilter PbiFilter::Intersection(std::vector<PbiFilter>&& filters)
{
    auto result = PbiFilter{ PbiFilter::INTERSECT };
    result.Add(std::move(filters));
    return result;
}

PbiFilter PbiFilter::Union(const std::vector<PbiFilter>& filters)
{
    auto result = PbiFilter{ PbiFilter::UNION };
    result.Add(filters);
    return result;
}

PbiFilter PbiFilter::Union(std::vector<PbiFilter>&& filters)
{
    auto result = PbiFilter{ PbiFilter::UNION };
    result.Add(std::move(filters));
    return result;
}
