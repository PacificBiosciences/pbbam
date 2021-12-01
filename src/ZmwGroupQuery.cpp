#include "PbbamInternalConfig.h"

#include <pbbam/ZmwGroupQuery.h>

#include <cstdint>

#include <algorithm>
#include <deque>

#include <pbbam/BamRecord.h>
#include <pbbam/CompositeBamReader.h>
#include <pbbam/PbiFilterTypes.h>

#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

class ZmwGroupQuery::ZmwGroupQueryPrivate
{
public:
    virtual ~ZmwGroupQueryPrivate() = default;
    virtual bool GetNext(std::vector<BamRecord>& records) = 0;

protected:
    ZmwGroupQueryPrivate() {}
};

class ZmwGroupQuery::WhitelistedZmwGroupQuery : public ZmwGroupQuery::ZmwGroupQueryPrivate
{
    using ReaderType = PbiFilterCompositeBamReader<Compare::Zmw>;

public:
    WhitelistedZmwGroupQuery(const std::vector<int32_t>& zmwWhitelist, const DataSet& dataset)
        : ZmwGroupQuery::ZmwGroupQueryPrivate()
        , whitelist_(zmwWhitelist.cbegin(), zmwWhitelist.cend())
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
    std::deque<int32_t> whitelist_;
    std::unique_ptr<ReaderType> reader_;
};

class ZmwGroupQuery::RoundRobinZmwGroupQuery : public ZmwGroupQuery::ZmwGroupQueryPrivate
{
public:
    RoundRobinZmwGroupQuery(const DataSet& dataset, const PbiFilter& pbiFilter)
        : ZmwGroupQuery::ZmwGroupQueryPrivate()
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

class ZmwGroupQuery::SequentialZmwGroupQuery : public ZmwGroupQuery::ZmwGroupQueryPrivate
{
public:
    SequentialZmwGroupQuery(const DataSet& dataset, const PbiFilter& pbiFilter)
        : ZmwGroupQuery::ZmwGroupQueryPrivate()
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

// ZmwGroupQuery(const DataSet& dataset, const PbiFilter& filter) {}

ZmwGroupQuery::ZmwGroupQuery(const DataSet& dataset, const PbiFilter& filter)
    : internal::IGroupQuery(), d_{std::make_unique<SequentialZmwGroupQuery>(dataset, filter)}
{
}

ZmwGroupQuery::ZmwGroupQuery(const std::vector<int32_t>& zmwWhitelist, const DataSet& dataset)
    : internal::IGroupQuery(), d_{std::make_unique<WhitelistedZmwGroupQuery>(zmwWhitelist, dataset)}
{
}

ZmwGroupQuery::ZmwGroupQuery(ZmwGroupQuery&&) noexcept = default;

ZmwGroupQuery& ZmwGroupQuery::operator=(ZmwGroupQuery&&) noexcept = default;

ZmwGroupQuery::~ZmwGroupQuery() = default;

bool ZmwGroupQuery::GetNext(std::vector<BamRecord>& records) { return d_->GetNext(records); }

}  // namespace BAM
}  // namespace PacBio
