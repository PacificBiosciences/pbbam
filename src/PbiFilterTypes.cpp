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
/// \file PbiFilterTypes.cpp
/// \brief Implements the built-in PBI filters.
//
// Author: Derek Barnett

#include "pbbam/PbiFilterTypes.h"
#include "StringUtils.h"
#include <boost/algorithm/string.hpp>
#include <sstream>
#include <string>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

template<typename T>
IndexList readLengthHelper(const std::vector<T>& start,
                           const std::vector<T>& end,
                           const T& value,
                           const Compare::Type cmp)
{
    assert(start.size() == end.size());

    auto result = IndexList{ };
    const auto numElements = start.size();
    for (size_t i = 0; i < numElements; ++i) {
        const auto readLength = end[i] - start[i];
        bool keep = false;
        switch(cmp) {
            case Compare::EQUAL              : keep = (readLength == value); break;
            case Compare::NOT_EQUAL          : keep = (readLength != value); break;
            case Compare::LESS_THAN          : keep = (readLength < value); break;
            case Compare::LESS_THAN_EQUAL    : keep = (readLength <= value); break;
            case Compare::GREATER_THAN       : keep = (readLength > value); break;
            case Compare::GREATER_THAN_EQUAL : keep = (readLength >= value); break;
            default:
                assert(false);
                throw std::runtime_error(string{"read length filter encountered unknown Compare::Type: "} +
                                         Compare::TypeToName(cmp));
        }

        if (keep)
            result.push_back(i);
    }
    return result;
}

static
PbiFilter filterFromMovieName(const string& movieName, bool includeCcs)
{
    // we'll match on any rgIds from our candidate list
    auto filter = PbiFilter{ PbiFilter::UNION };
    filter.Add(
    {
        PbiReadGroupFilter{ MakeReadGroupId(movieName, "POLYMERASE") },
        PbiReadGroupFilter{ MakeReadGroupId(movieName, "HQREGION") },
        PbiReadGroupFilter{ MakeReadGroupId(movieName, "SUBREAD") },
        PbiReadGroupFilter{ MakeReadGroupId(movieName, "SCRAP") },
        PbiReadGroupFilter{ MakeReadGroupId(movieName, "UNKNOWN") }
    });
    if (includeCcs)
        filter.Add(PbiReadGroupFilter{ MakeReadGroupId(movieName, "CCS") });

    return filter;
}

static
PbiFilter filterFromQueryName(const string& queryName)
{
    // split full name into moviename, holenumber
    const auto nameParts = internal::Split(queryName, '/');
    if (nameParts.size() != 3) {
        auto msg = string{ "PbiQueryNameFilter error: requested QNAME (" } + queryName;
        msg += string{ ") is not a valid PacBio BAM QNAME. See spec for details"};
        throw std::runtime_error(msg);
    }

    // main filter: {union of candidate rgIds} && zmw [&& qStart && qEnd](non-CCS reads)
    auto filter = PbiFilter{ };
    filter.Add(PbiZmwFilter{ stoi(nameParts.at(1)) }); // hole number

    const auto movieName = nameParts.at(0);

    // CCS (only 1 possible candidate rgId)
    if (nameParts.at(2) == "ccs")
        filter.Add(PbiReadGroupFilter{ MakeReadGroupId(movieName, "CCS") });

    // all other read types
    else {
        // we'll match on any read type that matches our qname
        // (except for CCS since it has a different QNAME anyway)
        const auto rgIdFilter = filterFromMovieName(movieName, false);
        filter.Add(rgIdFilter);

        // add qStart/qEnd filters to our main filter
        const auto queryIntervalParts = internal::Split(nameParts.at(2), '_');
        if (queryIntervalParts.size() != 2) {
            auto msg = string{ "PbiQueryNameFilter error: requested QNAME (" } + queryName;
            msg += string{ ") is not a valid PacBio BAM QNAME. See spec for details"};
            throw std::runtime_error(msg);
        }
        filter.Add(PbiQueryStartFilter{ stoi(queryIntervalParts.at(0)) });
        filter.Add(PbiQueryEndFilter{ stoi(queryIntervalParts.at(1)) });
    }
    return filter;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

// PbiAlignedLengthFilter

bool PbiAlignedLengthFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto& mappedData = idx.MappedData();
    const auto& aEnd    = mappedData.aEnd_.at(row) ;
    const auto& aStart  = mappedData.aStart_.at(row);
    const auto aLength = aEnd - aStart;
    return CompareHelper(aLength);
}

// PbiIdentityFilter

bool PbiIdentityFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto& mappedData = idx.MappedData();
    const auto& nMM  = mappedData.nMM_.at(row);
    const auto& nIndels = mappedData.NumDeletedAndInsertedBasesAt(row);
    const auto& nDel = nIndels.first;
    const auto& nIns = nIndels.second;

    const auto& basicData = idx.BasicData();
    const auto& qStart = basicData.qStart_.at(row);
    const auto& qEnd   = basicData.qEnd_.at(row);

    const auto readLength = qEnd - qStart;
    const auto nonMatches = nMM + nDel + nIns;
    const float identity  = 1.0 - (static_cast<float>(nonMatches)/static_cast<float>(readLength));

    return CompareHelper(identity);
}

// PbiMovieNameFilter

PbiMovieNameFilter::PbiMovieNameFilter(const std::string& movieName)
    : compositeFilter_(internal::filterFromMovieName(movieName, true)) // include CCS
{ }

PbiMovieNameFilter::PbiMovieNameFilter(const std::vector<std::string>& whitelist)
    : compositeFilter_(PbiFilter::UNION)
{
    for (const auto& movieName : whitelist)
        compositeFilter_.Add(internal::filterFromMovieName(movieName, true)); // include CCS
}

PbiMovieNameFilter::PbiMovieNameFilter(std::vector<std::string>&& whitelist)
    : compositeFilter_(PbiFilter::UNION)
{
    for (auto&& movieName : whitelist)
        compositeFilter_.Add(internal::filterFromMovieName(movieName, true)); // include CCS
}

// PbiQueryLengthFilter

bool PbiQueryLengthFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto& basicData = idx.BasicData();
    const auto& qStart = basicData.qStart_.at(row);
    const auto& qEnd   = basicData.qEnd_.at(row);
    const auto readLength = qEnd - qStart;
    return CompareHelper(readLength);
}

// PbiQueryNameFilter

PbiQueryNameFilter::PbiQueryNameFilter(const std::string& qname)
    : compositeFilter_(internal::filterFromQueryName(qname))
{ }

PbiQueryNameFilter::PbiQueryNameFilter(const std::vector<std::string>& whitelist)
    : compositeFilter_(PbiFilter::UNION)
{
    try {
        for (const auto& qname : whitelist)
            compositeFilter_.Add(internal::filterFromQueryName(qname));
    }
    // simply re-throw our own exception
    catch (std::runtime_error&) {
        throw;
    }
    // we may hit other exceptions (e.g. in stoi()) - but we'll pin on a bit of extra data
    catch (std::exception& e) {
        auto msg = string{ "PbiQueryNameFilter encountered error: " } + e.what();
        throw std::runtime_error(msg);
    }
}

PbiQueryNameFilter::PbiQueryNameFilter(std::vector<std::string>&& whitelist)
    : compositeFilter_(PbiFilter::UNION)
{
    try {
        for (const auto& qname : whitelist)
            compositeFilter_.Add(internal::filterFromQueryName(qname));
    }
    // simply re-throw our own exception
    catch (std::runtime_error&) {
        throw;
    }
    // we may hit other exceptions (e.g. in stoi()) - but we'll pin on a bit of extra data
    catch (std::exception& e) {
        auto msg = string{ "PbiQueryNameFilter encountered error: " } + e.what();
        throw std::runtime_error(msg);
    }
}

// PbiReferenceNameFilter

PbiReferenceNameFilter::PbiReferenceNameFilter(const std::string& rname,
                                               const Compare::Type cmp)
    : initialized_(false)
    , rname_(rname)
    , cmp_(cmp)
{
    if (cmp != Compare::EQUAL && cmp != Compare::NOT_EQUAL) {
        auto msg = std::string{ "Compare type: " };
        msg += Compare::TypeToName(cmp);
        msg += " not supported for PbiReferenceNameFilter (use one of Compare::EQUAL or Compare::NOT_EQUAL).";
        throw std::runtime_error(msg);
    }
}

PbiReferenceNameFilter::PbiReferenceNameFilter(const std::vector<std::string>& whitelist)
    : initialized_(false)
    , rnameWhitelist_(whitelist)
    , cmp_(Compare::EQUAL)
{ }

PbiReferenceNameFilter::PbiReferenceNameFilter(std::vector<std::string>&& whitelist)
    : initialized_(false)
    , rnameWhitelist_(std::move(whitelist))
    , cmp_(Compare::EQUAL)
{ }

bool PbiReferenceNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    if (!initialized_)
        Initialize(idx);
    return subFilter_.Accepts(idx, row);
}

void PbiReferenceNameFilter::Initialize(const PbiRawData& idx) const
{
    const auto pbiFilename = idx.Filename();
    const auto bamFilename = pbiFilename.substr(0, pbiFilename.length() - 4);
    const auto bamFile = BamFile{ bamFilename };

    // single-value
    if (rnameWhitelist_ == boost::none) {
        const auto tId = bamFile.ReferenceId(rname_);
        subFilter_ = PbiReferenceIdFilter{ tId, cmp_ };
    }

    // multi-value whitelist
    else {
        subFilter_ = PbiFilter(PbiFilter::UNION);
        for (const auto& rname : rnameWhitelist_.get())
            subFilter_.Add(PbiReferenceIdFilter{ bamFile.ReferenceId(rname) });
    }
    initialized_ = true;
}

