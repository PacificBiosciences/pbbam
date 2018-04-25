// File Description
/// \file PbiFilterTypes.cpp
/// \brief Implements the built-in PBI filters.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/PbiFilterTypes.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "StringUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T>
IndexList readLengthHelper(const std::vector<T>& start, const std::vector<T>& end, const T& value,
                           const Compare::Type cmp)
{
    assert(start.size() == end.size());

    auto result = IndexList{};
    const auto numElements = start.size();
    for (size_t i = 0; i < numElements; ++i) {
        const auto readLength = end[i] - start[i];
        bool keep = false;
        switch (cmp) {
            case Compare::EQUAL:
                keep = (readLength == value);
                break;
            case Compare::NOT_EQUAL:
                keep = (readLength != value);
                break;
            case Compare::LESS_THAN:
                keep = (readLength < value);
                break;
            case Compare::LESS_THAN_EQUAL:
                keep = (readLength <= value);
                break;
            case Compare::GREATER_THAN:
                keep = (readLength > value);
                break;
            case Compare::GREATER_THAN_EQUAL:
                keep = (readLength >= value);
                break;
            default:
                assert(false);
                throw std::runtime_error(
                    std::string{"read length filter encountered unknown Compare::Type: "} +
                    Compare::TypeToName(cmp));
        }

        if (keep) result.push_back(i);
    }
    return result;
}

static PbiFilter filterFromMovieName(const std::string& movieName, bool includeCcs)
{
    // we'll match on any rgIds from our candidate list
    auto filter = PbiFilter{PbiFilter::UNION};
    filter.Add({PbiReadGroupFilter{MakeReadGroupId(movieName, "POLYMERASE")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "HQREGION")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "SUBREAD")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "SCRAP")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "UNKNOWN")}});
    if (includeCcs) filter.Add(PbiReadGroupFilter{MakeReadGroupId(movieName, "CCS")});

    return filter;
}

}  // namespace internal

// PbiAlignedLengthFilter

bool PbiAlignedLengthFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto& mappedData = idx.MappedData();
    const auto& aEnd = mappedData.aEnd_.at(row);
    const auto& aStart = mappedData.aStart_.at(row);
    const auto aLength = aEnd - aStart;
    return CompareHelper(aLength);
}

// PbiIdentityFilter

bool PbiIdentityFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto& mappedData = idx.MappedData();
    const auto& nMM = mappedData.nMM_.at(row);
    const auto& nIndels = mappedData.NumDeletedAndInsertedBasesAt(row);
    const auto& nDel = nIndels.first;
    const auto& nIns = nIndels.second;

    const auto& basicData = idx.BasicData();
    const auto& qStart = basicData.qStart_.at(row);
    const auto& qEnd = basicData.qEnd_.at(row);

    const auto readLength = qEnd - qStart;
    const auto nonMatches = nMM + nDel + nIns;
    const float identity = 1.0 - (static_cast<float>(nonMatches) / static_cast<float>(readLength));

    return CompareHelper(identity);
}

// PbiMovieNameFilter

PbiMovieNameFilter::PbiMovieNameFilter(const std::string& movieName)
    : compositeFilter_(internal::filterFromMovieName(movieName, true))  // include CCS
{
}

PbiMovieNameFilter::PbiMovieNameFilter(const std::vector<std::string>& whitelist)
    : compositeFilter_(PbiFilter::UNION)
{
    for (const auto& movieName : whitelist)
        compositeFilter_.Add(internal::filterFromMovieName(movieName, true));  // include CCS
}

PbiMovieNameFilter::PbiMovieNameFilter(std::vector<std::string>&& whitelist)
    : compositeFilter_(PbiFilter::UNION)
{
    for (auto&& movieName : whitelist)
        compositeFilter_.Add(internal::filterFromMovieName(movieName, true));  // include CCS
}

// PbiQueryLengthFilter

bool PbiQueryLengthFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto& basicData = idx.BasicData();
    const auto& qStart = basicData.qStart_.at(row);
    const auto& qEnd = basicData.qEnd_.at(row);
    const auto readLength = qEnd - qStart;
    return CompareHelper(readLength);
}

// PbiQueryNameFilter

struct PbiQueryNameFilter::PbiQueryNameFilterPrivate
{
public:
    using QueryInterval = std::pair<int32_t, int32_t>;
    using QueryIntervals = std::set<QueryInterval>;
    using ZmwLookup = std::unordered_map<int32_t, QueryIntervals>;
    using ZmwLookupPtr = std::shared_ptr<ZmwLookup>;  // may be shared by more than one rgId
    using RgIdLookup = std::unordered_map<int32_t, ZmwLookupPtr>;

public:
    PbiQueryNameFilterPrivate(const std::vector<std::string>& whitelist)
    {
        for (const auto& queryName : whitelist) {

            if (queryName.find("transcript/") == 0)
                HandleName(queryName, RecordType::TRANSCRIPT);
            else if (queryName.find("/ccs") != std::string::npos)
                HandleName(queryName, RecordType::CCS);
            else
                HandleName(queryName, RecordType::UNKNOWN);
        }
    }

    PbiQueryNameFilterPrivate(const std::unique_ptr<PbiQueryNameFilterPrivate>& other)
    {
        if (other) lookup_ = other->lookup_;
    }

    bool Accepts(const PbiRawData& idx, const size_t row) const
    {
        const auto& basicData = idx.BasicData();

        // see if row's RGID known
        const auto& rgId = basicData.rgId_.at(row);
        const auto rgFound = lookup_.find(rgId);
        if (rgFound == lookup_.end()) return false;

        // see if row's ZMW known
        const auto& zmwPtr = rgFound->second;
        const auto zmw = basicData.holeNumber_.at(row);
        const auto zmwFound = zmwPtr->find(zmw);
        if (zmwFound == zmwPtr->end()) return false;

        // see if row's QueryStart/QueryEnd known
        // CCS names already covered in lookup construction phase
        const auto& queryIntervals = zmwFound->second;
        const auto qStart = basicData.qStart_.at(row);
        const auto qEnd = basicData.qEnd_.at(row);
        const auto queryInterval = std::make_pair(qStart, qEnd);
        return queryIntervals.find(queryInterval) != queryIntervals.end();
    }

    std::vector<int32_t> CandidateRgIds(const std::string& movieName, const RecordType type)
    {
        if (type == RecordType::CCS)
            return {ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "CCS"))};

        if (type == RecordType::TRANSCRIPT)
            return {ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "TRANSCRIPT"))};

        // we can't know for sure from QNAME alone
        return {ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "POLYMERASE")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "HQREGION")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "SUBREAD")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "SCRAP")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "UNKNOWN")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "ZMW"))};
    }

    void HandleName(const std::string& queryName, const RecordType type)
    {
        // split name into main parts
        const auto nameParts = internal::Split(queryName, '/');

        // verify syntax
        if (IsCcsOrTranscript(type)) {
            if (nameParts.size() != 2) {
                const auto typeName = (type == RecordType::CCS) ? "CCS" : "transcript";
                const auto msg = "PbiQueryNameFilter error: requested QNAME (" + queryName +
                                 ") is not valid for PacBio " + typeName +
                                 " reads. See spec for details.";
                throw std::runtime_error{msg};
            }
        } else {
            if (nameParts.size() != 3) {
                const auto msg = "PbiQueryNameFilter error: requested QNAME (" + queryName +
                                 ") is not a valid PacBio BAM QNAME. See spec for details";
                throw std::runtime_error{msg};
            }
        }

        // generate candidate read group IDs from movie name & record type, then
        // add to lookup table
        const auto zmwPtr = UpdateRgLookup(CandidateRgIds(nameParts.at(0), type));

        // add qStart/qEnd interval to zmw lookup
        const auto zmw = std::stoi(nameParts.at(1));
        if (IsCcsOrTranscript(type))
            UpdateZmwQueryIntervals(zmwPtr.get(), zmw, -1, -1);
        else {
            const auto queryIntervalParts = Split(nameParts.at(2), '_');
            if (queryIntervalParts.size() != 2) {
                auto msg = std::string{"PbiQueryNameFilter error: requested QNAME ("} + queryName;
                msg += std::string{") is not a valid PacBio BAM QNAME. See spec for details"};
                throw std::runtime_error{msg};
            }
            UpdateZmwQueryIntervals(zmwPtr.get(), zmw, std::stoi(queryIntervalParts.at(0)),
                                    std::stoi(queryIntervalParts.at(1)));
        }
    }

    ZmwLookupPtr UpdateRgLookup(std::vector<int32_t>&& rgIds)
    {
        assert(!rgIds.empty());

        ZmwLookupPtr zmwPtr;

        const auto rgFound = lookup_.find(rgIds.front());
        if (rgFound == lookup_.end()) {
            zmwPtr = std::make_shared<ZmwLookup>();
            for (const auto& rg : rgIds) {
                assert(lookup_.find(rg) == lookup_.end());
                lookup_.emplace(rg, zmwPtr);
            }
        } else {
#ifndef NDEBUG
            for (const auto& rg : rgIds)
                assert(lookup_.find(rg) != lookup_.end());
#endif
            zmwPtr = rgFound->second;
        }
        return zmwPtr;
    }

    // add QS/QE pair to ZMW lookup
    void UpdateZmwQueryIntervals(ZmwLookup* const zmwPtr, const int32_t zmw,
                                 const int32_t queryStart, const int32_t queryEnd)
    {
        const auto zmwFound = zmwPtr->find(zmw);
        if (zmwFound == zmwPtr->end()) zmwPtr->emplace(zmw, QueryIntervals{});
        auto& queryIntervals = zmwPtr->at(zmw);
        queryIntervals.emplace(std::make_pair(queryStart, queryEnd));
    }

private:
    RgIdLookup lookup_;
};

PbiQueryNameFilter::PbiQueryNameFilter(const std::string& qname)
    : d_(new PbiQueryNameFilter::PbiQueryNameFilterPrivate(std::vector<std::string>{1, qname}))
{
}

PbiQueryNameFilter::PbiQueryNameFilter(const std::vector<std::string>& whitelist)
    : d_(new PbiQueryNameFilter::PbiQueryNameFilterPrivate(whitelist))
{
}

PbiQueryNameFilter::PbiQueryNameFilter(const PbiQueryNameFilter& other)
    : d_(new PbiQueryNameFilter::PbiQueryNameFilterPrivate(other.d_))
{
}

PbiQueryNameFilter::~PbiQueryNameFilter() {}

bool PbiQueryNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    return d_->Accepts(idx, row);
}

// PbiReferenceNameFilter

PbiReferenceNameFilter::PbiReferenceNameFilter(std::string rname, Compare::Type cmp)
    : rname_(std::move(rname)), cmp_(cmp)
{
    if (cmp != Compare::EQUAL && cmp != Compare::NOT_EQUAL) {
        auto msg = std::string{"Compare type: "};
        msg += Compare::TypeToName(cmp);
        msg +=
            " not supported for PbiReferenceNameFilter (use one of Compare::EQUAL or "
            "Compare::NOT_EQUAL).";
        throw std::runtime_error(msg);
    }
}

PbiReferenceNameFilter::PbiReferenceNameFilter(const std::vector<std::string>& whitelist)
    : rnameWhitelist_(whitelist), cmp_(Compare::EQUAL)
{
}

PbiReferenceNameFilter::PbiReferenceNameFilter(std::vector<std::string>&& whitelist)
    : rnameWhitelist_(std::move(whitelist)), cmp_(Compare::EQUAL)
{
}

bool PbiReferenceNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    if (!initialized_) Initialize(idx);
    return subFilter_.Accepts(idx, row);
}

void PbiReferenceNameFilter::Initialize(const PbiRawData& idx) const
{
    const auto pbiFilename = idx.Filename();
    const auto bamFilename = pbiFilename.substr(0, pbiFilename.length() - 4);
    const auto bamFile = BamFile{bamFilename};

    // single-value
    if (rnameWhitelist_ == boost::none) {
        const auto tId = bamFile.ReferenceId(rname_);
        subFilter_ = PbiReferenceIdFilter{tId, cmp_};
    }

    // multi-value whitelist
    else {
        subFilter_ = PbiFilter(PbiFilter::UNION);
        for (const auto& rname : rnameWhitelist_.get())
            subFilter_.Add(PbiReferenceIdFilter{bamFile.ReferenceId(rname)});
    }
    initialized_ = true;
}

}  // namespace BAM
}  // namespace PacBio
