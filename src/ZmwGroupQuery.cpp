// File Description
/// \file ZmwQuery.cpp
/// \brief Implements the ZmwQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ZmwGroupQuery.h"

#include <algorithm>
#include <cstdint>
#include <deque>

#include "MemoryUtils.h"
#include "pbbam/BamRecord.h"
#include "pbbam/CompositeBamReader.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiFilterTypes.h"

namespace PacBio {
namespace BAM {

struct ZmwGroupQuery::ZmwGroupQueryPrivate
{
    using ReaderType = PbiFilterCompositeBamReader<Compare::Zmw>;

    ZmwGroupQueryPrivate(const std::vector<int32_t>& zmwWhitelist, const DataSet& dataset)
        : whitelist_(zmwWhitelist.cbegin(), zmwWhitelist.cend())
    {
        std::sort(whitelist_.begin(), whitelist_.end());
        whitelist_.erase(std::unique(whitelist_.begin(), whitelist_.end()), whitelist_.end());

        if (!whitelist_.empty()) {
            reader_ = std::make_unique<ReaderType>(PbiZmwFilter{whitelist_.front()}, dataset);
            whitelist_.pop_front();
        }
    }

    bool GetNext(std::vector<BamRecord>& records)
    {
        records.clear();
        if (!reader_) return false;

        // get all records matching ZMW
        BamRecord r;
        while (reader_->GetNext(r))
            records.push_back(r);

        // set next ZMW (if any left)
        if (!whitelist_.empty()) {
            reader_->Filter(PbiZmwFilter{whitelist_.front()});
            whitelist_.pop_front();
        }

        // otherwise destroy reader, next iteration will return false
        else
            reader_.reset();

        return true;
    }

    std::deque<int32_t> whitelist_;
    std::unique_ptr<ReaderType> reader_;
};

ZmwGroupQuery::ZmwGroupQuery(const std::vector<int32_t>& zmwWhitelist, const DataSet& dataset)
    : internal::IGroupQuery(), d_{std::make_unique<ZmwGroupQueryPrivate>(zmwWhitelist, dataset)}
{
}

ZmwGroupQuery::~ZmwGroupQuery() {}

bool ZmwGroupQuery::GetNext(std::vector<BamRecord>& records) { return d_->GetNext(records); }

}  // namespace BAM
}  // namespace PacBio
