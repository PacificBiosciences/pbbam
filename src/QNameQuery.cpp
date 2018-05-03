// File Description
/// \file QNameQuery.cpp
/// \brief Implements the QNameQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/QNameQuery.h"

#include <cassert>

#include <boost/optional.hpp>

#include "pbbam/CompositeBamReader.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

struct QNameQuery::QNameQueryPrivate
{
public:
    QNameQueryPrivate(const DataSet& dataset)
        : reader_{std::make_unique<SequentialCompositeBamReader>(dataset)}, nextRecord_(boost::none)
    {
    }

    bool GetNext(std::vector<BamRecord>& records)
    {
        records.clear();

        std::string groupRecordName;

        if (nextRecord_.is_initialized()) {
            BamRecord r = nextRecord_.get();
            groupRecordName = r.FullName();
            records.push_back(std::move(r));
            nextRecord_ = boost::none;
        }

        BamRecord record;
        while (reader_->GetNext(record)) {
            if (records.empty()) {
                groupRecordName = record.FullName();
                records.push_back(record);
            } else {
                assert(!records.empty());
                if (record.FullName() == groupRecordName)
                    records.push_back(record);
                else {
                    nextRecord_ = record;
                    return true;
                }
            }
        }
        return !records.empty();
    }

public:
    std::unique_ptr<SequentialCompositeBamReader> reader_;
    boost::optional<BamRecord> nextRecord_;
};

QNameQuery::QNameQuery(const DataSet& dataset)
    : internal::IGroupQuery(), d_{std::make_unique<QNameQueryPrivate>(dataset)}
{
}

QNameQuery::~QNameQuery() {}

bool QNameQuery::GetNext(std::vector<BamRecord>& records) { return d_->GetNext(records); }

}  // namespace BAM
}  // namespace PacBio
