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

#include "pbbam/MakeUnique.h"
#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {
namespace {

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
                throw std::runtime_error{"read length filter encountered unknown Compare::Type: " +
                                         Compare::TypeToName(cmp)};
        }

        if (keep) result.push_back(i);
    }
    return result;
}

PbiFilter filterFromMovieName(const std::string& movieName, bool includeCcs)
{
    //
    // All transcript-type reads (movieName == "transcript") have the same
    // read group ID. Calculate once & and create filters from that ID.
    //
    if (movieName == "transcript") {
        static const auto transcriptRgId = MakeReadGroupId("transcript", "TRANSCRIPT");
        return PbiFilter{PbiReadGroupFilter{transcriptRgId}};
    }

    //
    // For all other movie names, we can't determine read type up front, so we'll match
    // on any rgIds from a candidate list.
    //
    auto filter = PbiFilter{PbiFilter::UNION};
    filter.Add({PbiReadGroupFilter{MakeReadGroupId(movieName, "POLYMERASE")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "HQREGION")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "SUBREAD")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "SCRAP")},
                PbiReadGroupFilter{MakeReadGroupId(movieName, "UNKNOWN")}});
    if (includeCcs) filter.Add(PbiReadGroupFilter{MakeReadGroupId(movieName, "CCS")});

    return filter;
}

}  // anonymous

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

    const float readLength = qEnd - qStart;
    const float nonMatches = nMM + nDel + nIns;
    const float identity = 1.0f - (nonMatches / readLength);

    return CompareHelper(identity);
}

// PbiMovieNameFilter

PbiMovieNameFilter::PbiMovieNameFilter(const std::string& movieName, const Compare::Type cmp)
    : compositeFilter_{filterFromMovieName(movieName, true)}  // include CCS
    , cmp_{cmp}
{
}

PbiMovieNameFilter::PbiMovieNameFilter(const std::vector<std::string>& whitelist,
                                       const Compare::Type cmp)
    : compositeFilter_{PbiFilter::UNION}, cmp_{cmp}
{
    for (const auto& movieName : whitelist)
        compositeFilter_.Add(filterFromMovieName(movieName, true));  // include CCS
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

    PbiQueryNameFilterPrivate(const std::vector<std::string>& whitelist,
                              const Compare::Type cmp = Compare::EQUAL)
        : cmp_{cmp}
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
        if (other) {
            lookup_ = other->lookup_;
            cmp_ = other->cmp_;
        }
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

        const bool found = queryIntervals.find(queryInterval) != queryIntervals.end();
        if (cmp_ == Compare::EQUAL)
            return found;
        else if (cmp_ == Compare::NOT_EQUAL)
            return !found;
        else
            throw std::runtime_error{"unsupported compare type on query name filter"};
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
        const auto nameParts = Split(queryName, '/');

        // verify syntax
        if (IsCcsOrTranscript(type)) {
            if (nameParts.size() != 2) {
                const auto typeName = (type == RecordType::CCS) ? "CCS" : "transcript";
                throw std::runtime_error{"PbiQueryNameFilter error: requested QNAME (" + queryName +
                                         ") is not valid for PacBio " + typeName +
                                         " reads. See spec for details."};
            }
        } else {
            if (nameParts.size() != 3) {
                throw std::runtime_error{"PbiQueryNameFilter error: requested QNAME (" + queryName +
                                         ") is not a valid PacBio BAM QNAME. See spec for details"};
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
                throw std::runtime_error{"PbiQueryNameFilter error: requested QNAME (" + queryName +
                                         ") is not a valid PacBio BAM QNAME. See spec for details"};
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
    Compare::Type cmp_;
};

PbiQueryNameFilter::PbiQueryNameFilter(const std::string& qname, const Compare::Type cmp)
    : d_{std::make_unique<PbiQueryNameFilter::PbiQueryNameFilterPrivate>(
          std::vector<std::string>{1, qname}, cmp)}
{
}

PbiQueryNameFilter::PbiQueryNameFilter(const std::vector<std::string>& whitelist,
                                       const Compare::Type cmp)
    : d_{std::make_unique<PbiQueryNameFilter::PbiQueryNameFilterPrivate>(whitelist, cmp)}
{
}

PbiQueryNameFilter::PbiQueryNameFilter(const PbiQueryNameFilter& other)
    : d_{std::make_unique<PbiQueryNameFilter::PbiQueryNameFilterPrivate>(other.d_)}
{
}

PbiQueryNameFilter::~PbiQueryNameFilter() {}

bool PbiQueryNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    return d_->Accepts(idx, row);
}

// PbiReadGroupFilter

PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<int32_t>& whitelist,
                                       const Compare::Type cmp)
    : cmp_{cmp}
{
    if (cmp_ != Compare::EQUAL && cmp_ != Compare::NOT_EQUAL) {
        throw std::runtime_error{
            "Unsupported compare type for this property. "
            "Read group filter can only compare EQUAL or NOT_EQUAL"};
    }

    // Add RG ID & empty filter if not present. The empty filter will work for
    // non-barcoded IDs that match the expected number(s).
    //
    for (const auto& rgId : whitelist) {
        const auto found = lookup_.find(rgId);
        if (found == lookup_.cend()) lookup_.emplace(rgId, boost::none);
    }
}

PbiReadGroupFilter::PbiReadGroupFilter(const int32_t rgId, const Compare::Type cmp)
    : PbiReadGroupFilter{std::vector<int32_t>{rgId}, cmp}
{
}

PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<ReadGroupInfo>& whitelist,
                                       const Compare::Type cmp)
    : cmp_{cmp}
{
    if (cmp_ != Compare::EQUAL && cmp_ != Compare::NOT_EQUAL) {
        throw std::runtime_error{
            "Unsupported compare type for this property. "
            "Read group filter can only compare EQUAL or NOT_EQUAL"};
    }

    for (const auto& rg : whitelist) {
        // Add RG base ID with no filter if not present. The empty filter will
        // work for non-barcoded IDs. We'll add to it if the base read group ID
        // also has barcode labels,so that any barcode pair whitelisted for this
        // read group filter will be a match.
        //
        const auto idNum = ReadGroupInfo::IdToInt(rg.BaseId());
        const auto found = lookup_.find(idNum);
        if (found == lookup_.cend()) lookup_.emplace(idNum, boost::none);

        // Maybe add barcodes to base ID
        const auto barcodes = rg.Barcodes();
        if (barcodes) {
            const auto bcFor = static_cast<int16_t>(barcodes->first);
            const auto bcRev = static_cast<int16_t>(barcodes->second);
            auto& idBarcodes = lookup_.at(idNum);
            if (!idBarcodes) idBarcodes = std::vector<std::pair<int16_t, int16_t>>{};
            idBarcodes->push_back(std::make_pair(bcFor, bcRev));
        }
    }
}

PbiReadGroupFilter::PbiReadGroupFilter(const ReadGroupInfo& rg, const Compare::Type cmp)
    : PbiReadGroupFilter{std::vector<ReadGroupInfo>{rg}, cmp}
{
}

PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<std::string>& whitelist,
                                       const Compare::Type cmp)
{
    std::vector<ReadGroupInfo> readGroups;
    for (const auto rgId : whitelist)
        readGroups.push_back(rgId);
    *this = PbiReadGroupFilter{readGroups, cmp};
}

PbiReadGroupFilter::PbiReadGroupFilter(const std::string& rgId, const Compare::Type cmp)
    : PbiReadGroupFilter{ReadGroupInfo{rgId}, cmp}
{
}

bool PbiReadGroupFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto accepted = [this](const PbiRawData& index, const size_t i) {
        // Check that read group base ID is found.
        const auto rowRgId = index.BasicData().rgId_.at(i);
        const auto foundAt = lookup_.find(rowRgId);
        if (foundAt == lookup_.cend()) return false;

        // Read group's base ID is found, check for filtered barcodes.
        //
        // For non-barcoded read groups, the filter is empty. This is
        // essentially a no-op for allowing all candidate rows.
        //
        const auto& barcodes = foundAt->second;
        if (!barcodes) return true;

        // Return success on first match, otherwise no match found.
        for (const auto bcPair : *barcodes) {
            if (index.BarcodeData().bcForward_.at(i) == bcPair.first &&
                index.BarcodeData().bcReverse_.at(i) == bcPair.second) {
                return true;
            }
        }
        return false;

    }(idx, row);

    return (cmp_ == Compare::EQUAL ? accepted : !accepted);
}

// PbiReferenceNameFilter

PbiReferenceNameFilter::PbiReferenceNameFilter(std::string rname, Compare::Type cmp)
    : rname_{std::move(rname)}, cmp_{cmp}
{
    if (cmp != Compare::EQUAL && cmp != Compare::NOT_EQUAL) {
        throw std::runtime_error{
            "Compare type: " + Compare::TypeToName(cmp) +
            " not supported for PbiReferenceNameFilter (use one of Compare::EQUAL or "
            "Compare::NOT_EQUAL)."};
    }
}

PbiReferenceNameFilter::PbiReferenceNameFilter(std::vector<std::string> whitelist,
                                               const Compare::Type cmp)
    : rnameWhitelist_{std::move(whitelist)}, cmp_{cmp}
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
    const BamFile bamFile{bamFilename};

    // single-value
    if (rnameWhitelist_ == boost::none) {
        const auto tId = bamFile.ReferenceId(rname_);
        subFilter_ = PbiReferenceIdFilter{tId, cmp_};
    }

    // multi-value whitelist
    else {
        std::vector<int32_t> ids;
        for (const auto& rname : rnameWhitelist_.get())
            ids.push_back(bamFile.ReferenceId(rname));
        subFilter_ = PbiReferenceIdFilter{std::move(ids), cmp_};
    }
    initialized_ = true;
}

}  // namespace BAM
}  // namespace PacBio
