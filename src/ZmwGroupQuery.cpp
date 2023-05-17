#include "PbbamInternalConfig.h"

#include <pbbam/ZmwGroupQuery.h>

#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/CompositeBamReader.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiFilterTypes.h>
#include <pbbam/internal/QueryBase.h>
#include "MemoryUtils.h"

#include <algorithm>
#include <deque>
#include <memory>

#include <cstdint>

namespace PacBio {
namespace BAM {

class WhitelistedQuery : public internal::IGroupQuery
{
    using ReaderType = PbiFilterCompositeBamReader<Compare::Zmw>;

public:
    WhitelistedQuery(std::vector<std::int32_t> zmwWhitelist, const DataSet& dataset)
        : reader_{std::make_unique<ReaderType>(PbiZmwFilter{std::move(zmwWhitelist)}, dataset)}
    {
        if (!reader_->GetNext(currentRecord_)) {
            // No data is not an error, an actual error would be thrown by reader.
            // Reset reader to halt further iteration.
            reader_.reset();
        }
    }

    bool GetNext(std::vector<BamRecord>& records) override
    {
        records.clear();
        if (!reader_) {
            return false;
        }

        const int currentHoleNumber = currentRecord_.HoleNumber();
        records.push_back(currentRecord_);
        while (reader_->GetNext(currentRecord_)) {
            if (currentRecord_.HoleNumber() != currentHoleNumber) {
                // end of block, saving current record for next iteration
                return true;
            }
            // still in ZMW block
            records.push_back(currentRecord_);
        }

        // end of data, reset reader to halt further iteration
        reader_.reset();
        return true;
    }

private:
    BamRecord currentRecord_;
    std::unique_ptr<ReaderType> reader_;
};

///
/// Special case: aligned BAMs are not ordered by ZMW, but by mapping (chrom/pos).
/// This allows you to grab a block of ZMW subreads scattered anywhere in the file.
///
/// NOTE: This used to be the default behavior but is horribly inefficient for
///       blocks sorted by ZMW hole number.
///
class WhitelistedAlignmentQuery : public internal::IGroupQuery
{
    using ReaderType = PbiFilterCompositeBamReader<Compare::Zmw>;

public:
    WhitelistedAlignmentQuery(const std::vector<std::int32_t>& zmwWhitelist, const DataSet& dataset)
        : whitelist_(zmwWhitelist.cbegin(), zmwWhitelist.cend())
    {
        std::sort(whitelist_.begin(), whitelist_.end());
        whitelist_.erase(std::unique(whitelist_.begin(), whitelist_.end()), whitelist_.end());

        if (!whitelist_.empty()) {
            reader_ = std::make_unique<ReaderType>(PbiZmwFilter{whitelist_.front()}, dataset);
            whitelist_.pop_front();
        }
    }

    bool GetNext(std::vector<BamRecord>& records) override
    {
        records.clear();
        if (!reader_) {
            return false;
        }

        // get all records matching ZMW
        BamRecord r;
        while (reader_->GetNext(r)) {
            records.push_back(r);
        }

        // set next ZMW (if any left)
        if (!whitelist_.empty()) {
            reader_->Filter(PbiZmwFilter{whitelist_.front()});
            whitelist_.pop_front();
        }

        // otherwise destroy reader, next iteration will return false
        else {
            reader_.reset();
        }

        return true;
    }

private:
    std::deque<std::int32_t> whitelist_;
    std::unique_ptr<ReaderType> reader_;
};

std::unique_ptr<internal::IGroupQuery> MakeWhitelistedQuery(std::vector<std::int32_t> zmwWhitelist,
                                                            const DataSet& dataset)
{
    const auto mergedHeader = dataset.MergedHeader();
    if (mergedHeader.SortOrder() == "coordinate") {
        return std::make_unique<WhitelistedAlignmentQuery>(zmwWhitelist, dataset);
    }
    return std::make_unique<WhitelistedQuery>(std::move(zmwWhitelist), dataset);
}

class RoundRobinZmwGroupQuery : public internal::IGroupQuery
{
public:
    RoundRobinZmwGroupQuery(const DataSet& dataset, const PbiFilter& pbiFilter)
    {
        const auto bamFilenames = dataset.BamFilenames();
        for (const auto& fn : bamFilenames) {
            // create reader for file
            auto makeReader = [&]() -> std::unique_ptr<BamReader> {
                if (pbiFilter.IsEmpty()) {
                    return std::make_unique<BamReader>(fn);
                } else {
                    return std::make_unique<PbiIndexedBamReader>(pbiFilter, fn);
                }
            };
            internal::CompositeMergeItem item{makeReader()};

            // try load first record, ignore file if nothing found
            if (item.reader->GetNext(item.record)) {
                readerItems_.push_back(std::move(item));
            }
        }
    }

    bool GetNext(std::vector<BamRecord>& records) override
    {
        records.clear();

        // quick exit if nothing left
        if (readerItems_.empty()) {
            return false;
        }

        // pop first reader from the queue & store its record
        auto firstIter = readerItems_.begin();
        internal::CompositeMergeItem item{std::move(firstIter->reader),
                                          std::move(firstIter->record)};
        readerItems_.pop_front();
        auto zmw = item.record.HoleNumber();
        records.push_back(item.record);

        // try to read all records from the current reader, matching ZMW numbers
        while (true) {
            if (item.reader->GetNext(item.record)) {

                // if same ZMW, store and continue
                // else stop reading & re-queue this reader for later
                if (item.record.HoleNumber() == zmw) {
                    records.push_back(item.record);
                } else {
                    readerItems_.push_back(std::move(item));
                    break;
                }
            }

            // no data remaining for this reader, let it go
            else {
                break;
            }
        }

        return !records.empty();
    }

private:
    std::deque<internal::CompositeMergeItem> readerItems_;
};

class SequentialZmwGroupQuery : public internal::IGroupQuery
{
public:
    SequentialZmwGroupQuery(const DataSet& dataset, const PbiFilter& pbiFilter)
    {
        const auto bamFilenames = dataset.BamFilenames();
        for (const auto& fn : bamFilenames) {
            // create reader for file
            auto makeReader = [&]() -> std::unique_ptr<BamReader> {
                if (pbiFilter.IsEmpty()) {
                    return std::make_unique<BamReader>(fn);
                } else {
                    return std::make_unique<PbiIndexedBamReader>(pbiFilter, fn);
                }
            };
            internal::CompositeMergeItem item{makeReader()};

            // try load first record, ignore file if nothing found
            if (item.reader->GetNext(item.record)) {
                readerItems_.push_back(std::move(item));
            }
        }
    }

    bool GetNext(std::vector<BamRecord>& records) override
    {
        records.clear();

        // quick exit if nothing left
        if (readerItems_.empty()) {
            return false;
        }

        // pop first reader from the queue & store its record
        auto firstIter = readerItems_.begin();
        internal::CompositeMergeItem item{std::move(firstIter->reader),
                                          std::move(firstIter->record)};
        readerItems_.pop_front();
        auto zmw = item.record.HoleNumber();
        records.push_back(item.record);

        while (true) {
            if (item.reader->GetNext(item.record)) {
                // if same ZMW, store and continue
                // else stop reading and keep this reader/record for next time
                if (item.record.HoleNumber() == zmw) {
                    records.push_back(item.record);
                } else {
                    readerItems_.push_front(std::move(item));
                    break;
                }
            }

            // no data remaining for this reader, let it go
            else {
                break;
            }
        }
        return !records.empty();
    }

private:
    std::deque<internal::CompositeMergeItem> readerItems_;
};

ZmwGroupQuery::ZmwGroupQuery(const DataSet& dataset, const ZmwFileIterationMode iterationMode,
                             const DataSetFilterMode filterMode)
    : internal::IGroupQuery()
{
    PbiFilter filter;
    if (filterMode == DataSetFilterMode::APPLY) {
        filter = PbiFilter::FromDataSet(dataset);
    }

    if (iterationMode == ZmwFileIterationMode::SEQUENTIAL) {
        d_ = std::make_unique<SequentialZmwGroupQuery>(dataset, filter);
    } else {
        d_ = std::make_unique<RoundRobinZmwGroupQuery>(dataset, filter);
    }
}

ZmwGroupQuery::ZmwGroupQuery(const DataSet& dataset, const PbiFilter& filter)
    : internal::IGroupQuery(), d_{std::make_unique<SequentialZmwGroupQuery>(dataset, filter)}
{}

ZmwGroupQuery::ZmwGroupQuery(std::vector<std::int32_t> zmwWhitelist, const DataSet& dataset)
    : internal::IGroupQuery(), d_{MakeWhitelistedQuery(std::move(zmwWhitelist), dataset)}
{}

ZmwGroupQuery::ZmwGroupQuery(ZmwGroupQuery&&) noexcept = default;

ZmwGroupQuery& ZmwGroupQuery::operator=(ZmwGroupQuery&&) noexcept = default;

ZmwGroupQuery::~ZmwGroupQuery() = default;

bool ZmwGroupQuery::GetNext(std::vector<BamRecord>& records) { return d_->GetNext(records); }

}  // namespace BAM
}  // namespace PacBio
