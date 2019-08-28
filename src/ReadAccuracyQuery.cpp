// File Description
/// \file ReadAccuracyQuery.cpp
/// \brief Implements the ReadAccuracyQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ReadAccuracyQuery.h"

#include "pbbam/CompositeBamReader.h"
#include "pbbam/PbiFilterTypes.h"

namespace PacBio {
namespace BAM {

class ReadAccuracyQuery::ReadAccuracyQueryPrivate
{
public:
    ReadAccuracyQueryPrivate(const Accuracy accuracy, const Compare::Type compareType,
                             const DataSet& dataset)
        : reader_{PbiReadAccuracyFilter{accuracy, compareType}, dataset}
    {
    }

    PbiFilterCompositeBamReader<Compare::None> reader_;  // unsorted
};

ReadAccuracyQuery::ReadAccuracyQuery(const Accuracy accuracy, const Compare::Type compareType,
                                     const DataSet& dataset)
    : internal::IQuery()
    , d_{std::make_unique<ReadAccuracyQueryPrivate>(accuracy, compareType, dataset)}
{
}

ReadAccuracyQuery::~ReadAccuracyQuery() = default;

bool ReadAccuracyQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

uint32_t ReadAccuracyQuery::NumReads() const { return d_->reader_.NumReads(); }

}  // namespace BAM
}  // namespace PacBio
