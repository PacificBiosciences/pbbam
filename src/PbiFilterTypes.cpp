#include "PbbamInternalConfig.h"

#include <pbbam/PbiFilterTypes.h>

#include <pbbam/StringUtilities.h>

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <optional>
#include <string>
#include <unordered_map>

#include <cassert>
#include <cstddef>
#include <cstdint>

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
                throw std::runtime_error{
                    "[pbbam] PBI filter ERROR: read length filter encountered unknown compare "
                    "type: " +
                    Compare::TypeToName(cmp)};
        }

        if (keep) {
            result.push_back(i);
        }
    }
    return result;
}

}  // namespace

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
    : PbiMovieNameFilter{{1, movieName}, cmp}
{
}

PbiMovieNameFilter::PbiMovieNameFilter(const std::vector<std::string>& movieNames,
                                       const Compare::Type cmp)
    : cmp_{cmp}
{
    if (cmp_ == Compare::EQUAL) {
        cmp_ = Compare::CONTAINS;
    } else if (cmp_ == Compare::NOT_EQUAL) {
        cmp_ = Compare::NOT_CONTAINS;
    }

    if (cmp_ != Compare::CONTAINS && cmp_ != Compare::NOT_CONTAINS) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: unsupported compare type (" + Compare::TypeToName(cmp) +
            ") for this property. "
            "Movie name filter can only compare equality or presence in whitelist/blacklist."};
    }

    // clang-format off
    for (const auto& movieName : movieNames) {
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "CCS")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "TRANSCRIPT")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "POLYMERASE")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "HQREGION")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "SUBREAD")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "SCRAP")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "UNKNOWN")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "ZMW")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "CCS")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "TRANSCRIPT")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "POLYMERASE")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "HQREGION")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "SUBREAD")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "SCRAP")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "UNKNOWN")));
        candidateRgIds_.insert(ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "ZMW")));
        movieNames_.insert(movieName);
    }
    // clang-format on
}

bool PbiMovieNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto accepted = [this](const PbiRawData& index, const size_t i) {

        // straightforward lookup
        const auto& rgId = index.BasicData().rgId_.at(i);
        const auto foundAt = candidateRgIds_.find(rgId);
        if (foundAt != candidateRgIds_.cend()) {
            return true;
        }

        // if no barcode context available, record movie name fails
        if (!index.HasBarcodeData()) {
            return false;
        }

        // try barcoded RG IDs
        const auto& barcodeData = index.BarcodeData();
        const auto barcodes =
            std::make_pair(barcodeData.bcForward_.at(i), barcodeData.bcReverse_.at(i));
        for (const auto& movieName : movieNames_) {
            const auto tryBarcodedType = [&](const std::string& readType) {
                int32_t barcodedId =
                    ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, readType, barcodes));
                if (barcodedId == rgId) {
                    candidateRgIds_.insert(barcodedId);  // found combo, save for future lookup
                    return true;
                }

                barcodedId =
                    ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, readType, barcodes));
                if (barcodedId == rgId) {
                    candidateRgIds_.insert(barcodedId);  // found combo, save for future lookup
                    return true;
                }

                return false;
            };

            if (tryBarcodedType("CCS")) {
                return true;
            }
            if (tryBarcodedType("TRANSCRIPT")) {
                return true;
            }
            if (tryBarcodedType("SUBREAD")) {
                return true;
            }
            if (tryBarcodedType("ZMW")) {
                return true;
            }
            if (tryBarcodedType("POLYMERASE")) {
                return true;
            }
            if (tryBarcodedType("HQREGION")) {
                return true;
            }
            if (tryBarcodedType("SCRAP")) {
                return true;
            }
            if (tryBarcodedType("UNKNOWN")) {
                return true;
            }
        }

        // not found at all
        return false;
    }(idx, row);

    assert(cmp_ == Compare::CONTAINS || cmp_ == Compare::NOT_CONTAINS);
    return (cmp_ == Compare::CONTAINS ? accepted : !accepted);
}

// PbiNumSubreadsFilter

struct PbiNumSubreadsFilter::PbiNumSubreadsFilterPrivate
{
    PbiNumSubreadsFilterPrivate(int numSubreads, const Compare::Type cmp)
        : numSubreads_{numSubreads}, cmp_{cmp}
    {
    }

    PbiNumSubreadsFilterPrivate(const std::unique_ptr<PbiNumSubreadsFilterPrivate>& other)
    {
        if (other) {
            numSubreads_ = other->numSubreads_;
            lookup_ = other->lookup_;
            cmp_ = other->cmp_;
        }
    }

    bool Accepts(const PbiRawData& idx, const size_t row) const
    {
        // lazy-load
        if (!lookup_) {
            InitializeLookup(idx);
        }
        const auto holeNumber = idx.BasicData().holeNumber_[row];
        return (lookup_->find(holeNumber) != lookup_->cend());
    }

    void InitializeLookup(const PbiRawData& idx) const
    {
        lookup_ = std::set<int32_t>{};
        const auto& zmws = idx.BasicData().holeNumber_;

        auto shouldKeep = [this](int count) {
            return Compare::Check(count, this->numSubreads_, this->cmp_);
        };

        auto start = zmws.begin();
        auto current = zmws.begin();
        int count = 0;
        while (current < zmws.cend()) {
            if (*start != *current) {
                if (shouldKeep(count)) {
                    lookup_->insert(*start);
                }
                start = current;
                count = 0;
            }
            ++count;
            ++current;
        }
        if (shouldKeep(count)) {
            lookup_->insert(*start);
        }
    }

    int numSubreads_;
    Compare::Type cmp_;
    mutable std::optional<std::set<int32_t>> lookup_;  // mutable for lazy-load
};

PbiNumSubreadsFilter::PbiNumSubreadsFilter(int numSubreads, const Compare::Type cmp)
    : d_{std::make_unique<PbiNumSubreadsFilter::PbiNumSubreadsFilterPrivate>(numSubreads, cmp)}
{
}

PbiNumSubreadsFilter::PbiNumSubreadsFilter(const PbiNumSubreadsFilter& other)
    : d_{std::make_unique<PbiNumSubreadsFilter::PbiNumSubreadsFilterPrivate>(other.d_)}
{
}

PbiNumSubreadsFilter::PbiNumSubreadsFilter(PbiNumSubreadsFilter&&) noexcept = default;

PbiNumSubreadsFilter& PbiNumSubreadsFilter::operator=(PbiNumSubreadsFilter&&) noexcept = default;

PbiNumSubreadsFilter::~PbiNumSubreadsFilter() = default;

bool PbiNumSubreadsFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    return d_->Accepts(idx, row);
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
    using ZmwData = std::unordered_map<int32_t, std::optional<QueryIntervals>>;
    using RgIdLookup = std::unordered_map<int32_t, std::shared_ptr<ZmwData>>;

    PbiQueryNameFilterPrivate(const std::vector<std::string>& queryNames,
                              const Compare::Type cmp = Compare::EQUAL)
        : cmp_{cmp}
    {
        for (const auto& queryName : queryNames) {

            if (queryName.find("transcript/") == 0) {
                HandleName(queryName, RecordType::TRANSCRIPT);
            } else if (queryName.find("/ccs") != std::string::npos) {
                HandleName(queryName, RecordType::CCS);
            } else {
                HandleName(queryName, RecordType::UNKNOWN);
            }
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

        const bool found = [&]() {
            // see if row's RGID known
            const auto& rgId = basicData.rgId_.at(row);
            const auto rgFound = lookup_.find(rgId);
            if (rgFound == lookup_.end()) {
                return false;
            }

            // see if row's ZMW known
            const auto& zmwPtr = rgFound->second;
            const auto zmw = basicData.holeNumber_.at(row);
            const auto zmwFound = zmwPtr->find(zmw);
            if (zmwFound == zmwPtr->end()) {
                return false;
            }

            // see if row's QueryStart/QueryEnd known
            // CCS names already covered in lookup construction phase
            const auto& queryIntervals = zmwFound->second;
            if (queryIntervals) {
                const auto qStart = basicData.qStart_.at(row);
                const auto qEnd = basicData.qEnd_.at(row);
                const QueryInterval queryInterval = std::make_pair(qStart, qEnd);
                return (queryIntervals->find(queryInterval) != queryIntervals->end());
            } else {
                // CCS or transcript record
                return true;
            }
        }();

        if (cmp_ == Compare::EQUAL || cmp_ == Compare::CONTAINS) {
            return found;
        } else if (cmp_ == Compare::NOT_EQUAL || cmp_ == Compare::NOT_CONTAINS) {
            return !found;
        } else {
            throw std::runtime_error{
                "[pbbam] PBI filter ERROR: unsupported compare type on query name filter"};
        }
    }

    std::vector<int32_t> CandidateRgIds(const std::string& movieName, const RecordType type)
    {
        if (type == RecordType::CCS) {
            return {
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "CCS")),
                ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "CCS")),
            };
        }

        if (type == RecordType::TRANSCRIPT) {
            return {ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "TRANSCRIPT")),
                    ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "TRANSCRIPT"))};
        }

        // we can't know for sure from QNAME alone
        return {ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "POLYMERASE")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "HQREGION")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "SUBREAD")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "SCRAP")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "UNKNOWN")),
                ReadGroupInfo::IdToInt(MakeReadGroupId(movieName, "ZMW")),
                ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "POLYMERASE")),
                ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "HQREGION")),
                ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "SUBREAD")),
                ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "SCRAP")),
                ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "UNKNOWN")),
                ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(movieName, "ZMW"))};
    }

    void HandleName(const std::string& queryName, const RecordType type)
    {
        const auto nameParts = Split(queryName, '/');
        if (nameParts.size() < 2) {
            throw std::runtime_error{"[pbbam] PBI filter ERROR: requested QNAME (" + queryName +
                                     ") is not a valid PacBio BAM QNAME. See spec for details"};
        }

        // generate candidate read group IDs from movie name & record type, then
        // add to lookup table
        const std::shared_ptr<ZmwData> zmw = UpdateRgLookup(CandidateRgIds(nameParts.at(0), type));

        // add ZMW to read group. Add qStart/qEnd to ZMW if not a CCS/transcript record
        const auto zmwId = [&]() {
            try {
                return std::stoi(nameParts.at(1));
            } catch (const std::invalid_argument&) {
                throw std::runtime_error{
                    "[pbbam] PBI filter ERROR: requested QNAME (" + queryName +
                    ") is not a valid PacBio BAM QNAME. ZMW id must be a number."};
            }
        }();

        if (IsCcsOrTranscript(type)) {
            zmw->emplace(zmwId, std::optional<QueryIntervals>{});
        } else {

            const auto queryIntervalParts = Split(nameParts.at(2), '_');
            if (queryIntervalParts.size() != 2) {
                throw std::runtime_error{"[pbbam] PBI filter ERROR: requested QNAME (" + queryName +
                                         ") is not a valid PacBio BAM QNAME. See spec for details"};
            }

            const auto queryInterval = [&]() {
                try {
                    return std::make_pair(std::stoi(queryIntervalParts.at(0)),
                                          std::stoi(queryIntervalParts.at(1)));
                } catch (const std::invalid_argument&) {
                    throw std::runtime_error{
                        "[pbbam] PBI filter ERROR: requested QNAME (" + queryName +
                        ") is not a valid PacBio BAM QNAME. qStart/qEnd must be numbers."};
                }
            }();

            const auto zmwResult = zmw->emplace(zmwId, QueryIntervals{});
            const auto zmwIter = zmwResult.first;
            auto& queryIntervals = zmwIter->second;
            queryIntervals->emplace(queryInterval);
        }
    }

    std::shared_ptr<ZmwData> UpdateRgLookup(const std::vector<int32_t>& rgIds)
    {
        assert(!rgIds.empty());

        std::shared_ptr<ZmwData> zmw;

        const auto rgFound = lookup_.find(rgIds.front());
        if (rgFound == lookup_.end()) {
            zmw = std::make_shared<ZmwData>();
            for (const auto& rg : rgIds) {
                // Extra RG hashes (fixed & legacy) have been calculated as
                // candidates, but sometimes these are the same. Only store once.
                if (lookup_.find(rg) == lookup_.end()) {
                    lookup_.emplace(rg, zmw);
                }
            }
        } else {
#ifndef NDEBUG
            for (const auto& rg : rgIds) {
                assert(lookup_.find(rg) != lookup_.end());
            }
#endif
            zmw = rgFound->second;
        }
        return zmw;
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

PbiQueryNameFilter::PbiQueryNameFilter(const std::vector<std::string>& queryNames,
                                       const Compare::Type cmp)
    : d_{std::make_unique<PbiQueryNameFilter::PbiQueryNameFilterPrivate>(queryNames, cmp)}
{
}

PbiQueryNameFilter::PbiQueryNameFilter(const PbiQueryNameFilter& other)
    : d_{std::make_unique<PbiQueryNameFilter::PbiQueryNameFilterPrivate>(other.d_)}
{
}

PbiQueryNameFilter::PbiQueryNameFilter(PbiQueryNameFilter&&) noexcept = default;

PbiQueryNameFilter& PbiQueryNameFilter::operator=(PbiQueryNameFilter&&) noexcept = default;

PbiQueryNameFilter::~PbiQueryNameFilter() = default;

bool PbiQueryNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    return d_->Accepts(idx, row);
}

// PbiReadGroupFilter

PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<int32_t>& rgIds, const Compare::Type cmp)
    : cmp_{cmp}
{
    std::vector<ReadGroupInfo> readGroups{rgIds.size()};
    std::transform(rgIds.cbegin(), rgIds.cend(), readGroups.begin(),
                   [](const int32_t rgId) { return ReadGroupInfo{ReadGroupInfo::IntToId(rgId)}; });
    AddReadGroups(readGroups);
}

PbiReadGroupFilter::PbiReadGroupFilter(const int32_t rgId, const Compare::Type cmp)
    : PbiReadGroupFilter{std::vector<int32_t>{rgId}, cmp}
{
}

PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<ReadGroupInfo>& readGroups,
                                       const Compare::Type cmp)
    : cmp_{cmp}
{
    AddReadGroups(readGroups);
}

PbiReadGroupFilter::PbiReadGroupFilter(const ReadGroupInfo& rg, const Compare::Type cmp)
    : PbiReadGroupFilter{std::vector<ReadGroupInfo>{rg}, cmp}
{
}

PbiReadGroupFilter::PbiReadGroupFilter(const std::vector<std::string>& rgIds,
                                       const Compare::Type cmp)
    : cmp_{cmp}
{
    std::vector<ReadGroupInfo> readGroups{rgIds.size()};
    std::transform(rgIds.cbegin(), rgIds.cend(), readGroups.begin(),
                   [](const std::string& rgId) { return ReadGroupInfo{rgId}; });
    AddReadGroups(readGroups);
}

PbiReadGroupFilter::PbiReadGroupFilter(const std::string& rgId,
                                       const Compare::Type cmp)  //: cmp_{cmp}
    : PbiReadGroupFilter{std::vector<std::string>{rgId}, cmp}
{
}

bool PbiReadGroupFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto DoFiltersMatch = [&](const int32_t rowRgId) {
        const auto foundInFilterList = readGroups_.find(rowRgId);
        if (foundInFilterList == readGroups_.cend()) {
            return false;
        }

        // matching ID found, check for potential barcode requirements

        if (idx.HasBarcodeData()) {
            const int16_t rowBcForward = idx.BarcodeData().bcForward_.at(row);
            const int16_t rowBcReverse = idx.BarcodeData().bcReverse_.at(row);

            for (const auto& filterReadGroup : foundInFilterList->second) {
                const auto& filterBarcodes = filterReadGroup.Barcodes();
                if (filterBarcodes) {
                    const int16_t filterBcForward = filterBarcodes->first;
                    const int16_t filterBcReverse = filterBarcodes->second;
                    if ((rowBcForward == filterBcForward) && (rowBcReverse == filterBcReverse)) {
                        // found matching barcodes
                        return true;
                    }
                }
            }

            // no read groups in filter match this index row's barcodes
            return false;
        } else {
            for (const auto& filterReadGroup : foundInFilterList->second) {
                const auto& filterBarcodes = filterReadGroup.Barcodes();
                if (!filterBarcodes) {
                    // found a read group that matches ID & does not require a barcode match
                    return true;
                }
            }

            // all filter read groups require barcodes, but index does not
            // contain any barcode information
            return false;
        }
    };

    const int rowRgId = idx.BasicData().rgId_.at(row);
    const bool rowMatched = DoFiltersMatch(rowRgId);
    const bool lookingForEquality = (cmp_ == Compare::CONTAINS) || (cmp_ == Compare::EQUAL);
    return (lookingForEquality ? rowMatched : !rowMatched);
}

void PbiReadGroupFilter::AddReadGroups(const std::vector<ReadGroupInfo>& readGroups)
{
    if (cmp_ == Compare::EQUAL) {
        cmp_ = Compare::CONTAINS;
    } else if (cmp_ == Compare::NOT_EQUAL) {
        cmp_ = Compare::NOT_CONTAINS;
    }

    if (cmp_ != Compare::CONTAINS && cmp_ != Compare::NOT_CONTAINS) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: unsupported compare type (" +
                                 Compare::TypeToName(cmp_) +
                                 ") for this property. Read group filter can only compare equality "
                                 "or presence in whitelist/blacklist."};
    }

    //
    // Ensure we track all potential representations of a read group's ID.
    //
    // NOTE: Storing the read group object more than once for equivalent IDs is
    //       allowed here. The matching phase does a linear walk over the read groups
    //       stored here to determine a match. This does not change the result.
    //
    for (const auto& rg : readGroups) {
        const std::string rgId = rg.Id();
        readGroups_[ReadGroupInfo::IdToInt(rgId)].push_back(rg);
        readGroups_[ReadGroupInfo::IdToInt(ReadGroupInfo::GetBaseId(rgId))].push_back(rg);
        readGroups_[ReadGroupInfo::IdToInt(MakeReadGroupId(rg))].push_back(rg);
        readGroups_[ReadGroupInfo::IdToInt(MakeLegacyReadGroupId(rg))].push_back(rg);
    }
}

// PbiReferenceNameFilter

PbiReferenceNameFilter::PbiReferenceNameFilter(std::string rname, Compare::Type cmp)
    : rname_{std::move(rname)}, cmp_{cmp}
{
    Validate();
}

PbiReferenceNameFilter::PbiReferenceNameFilter(std::vector<std::string> rnames,
                                               const Compare::Type cmp)
    : rnameWhitelist_{std::move(rnames)}, cmp_{cmp}
{
    Validate();
}

bool PbiReferenceNameFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    if (!initialized_) {
        Initialize(idx);
    }
    return subFilter_.Accepts(idx, row);
}

void PbiReferenceNameFilter::Initialize(const PbiRawData& idx) const
{

    // fetch BAM header info associate with this index
    const auto pbiFilename = idx.Filename();
    const auto bamFilename = pbiFilename.substr(0, pbiFilename.length() - 4);
    const BamFile bamFile{bamFilename};

    // single-value
    if (!rnameWhitelist_) {
        const auto tId = bamFile.ReferenceId(rname_);
        subFilter_ = PbiReferenceIdFilter{tId, cmp_};
    }

    // multi-value (whitelist/blacklist)
    else {
        std::vector<int32_t> ids;
        for (const auto& rname : *rnameWhitelist_) {
            ids.push_back(bamFile.ReferenceId(rname));
        }
        subFilter_ = PbiReferenceIdFilter{std::move(ids), cmp_};
    }
    initialized_ = true;
}

void PbiReferenceNameFilter::Validate() const
{
    // double-check valid compare type
    const bool compareTypeOk = [&]() {
        if (cmp_ == Compare::EQUAL) {
            return true;
        }
        if (cmp_ == Compare::NOT_EQUAL) {
            return true;
        }
        if (cmp_ == Compare::CONTAINS) {
            return true;
        }
        if (cmp_ == Compare::NOT_CONTAINS) {
            return true;
        }
        return false;
    }();
    if (!compareTypeOk) {
        throw std::runtime_error{"[pbbam] PBI filter ERROR: unsupported compare type (" +
                                 Compare::TypeToName(cmp_) +
                                 ") for this property. "
                                 "Reference name filter can only compare equality or presence "
                                 "in whitelist/blacklist."};
    }
}

// PbiZmwFilter

PbiZmwFilter::PbiZmwFilter(const int32_t zmw, const Compare::Type cmp) : cmp_{cmp}, singleZmw_{zmw}
{
}

PbiZmwFilter::PbiZmwFilter(std::vector<int32_t> whitelist, const Compare::Type cmp)
    : cmp_{cmp}, zmwLookup_{whitelist.begin(), whitelist.end()}
{
    // verify white list compare type
    if (cmp_ == Compare::EQUAL) {
        cmp_ = Compare::CONTAINS;
    } else if (cmp_ == Compare::NOT_EQUAL) {
        cmp_ = Compare::NOT_CONTAINS;
    }
    if (cmp_ != Compare::CONTAINS && cmp_ != Compare::NOT_CONTAINS) {
        throw std::runtime_error{
            "[pbbam] PBI filter ERROR: multi-valued filters (e.g. whitelists) can only check "
            "for "
            "containment."};
    }
}

bool PbiZmwFilter::Accepts(const PbiRawData& idx, const size_t row) const
{
    const auto zmw = idx.BasicData().holeNumber_.at(row);
    if (cmp_ == Compare::CONTAINS) {
        return zmwLookup_.find(zmw) != zmwLookup_.cend();
    } else if (cmp_ == Compare::NOT_CONTAINS) {
        return zmwLookup_.find(zmw) == zmwLookup_.cend();
    } else {
        return Compare::Check(zmw, singleZmw_, cmp_);
    }
}

}  // namespace BAM
}  // namespace PacBio
