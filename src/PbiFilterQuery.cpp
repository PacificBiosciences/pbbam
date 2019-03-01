// File Description
/// \file PbiFilterQuery.cpp
/// \brief Implements the PbiFilterQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/PbiFilterQuery.h"

#include <iostream>

#include "pbbam/CompositeBamReader.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

class PbiFilterQuery::PbiFilterQueryPrivate
{
public:
    PbiFilterQueryPrivate(const PbiFilter& filter, const DataSet& dataset,
                          const PbiIndexCache& cache)
        : reader_{filter, dataset, cache}
    {
    }

    PbiFilterCompositeBamReader<Compare::None> reader_;  // unsorted
};

PbiFilterQuery::PbiFilterQuery(const DataSet& dataset)
    : PbiFilterQuery{PbiFilter::FromDataSet(dataset), dataset, MakePbiIndexCache(dataset)}
{
}

PbiFilterQuery::PbiFilterQuery(const DataSet& dataset, const PbiIndexCache& cache)
    : PbiFilterQuery{PbiFilter::FromDataSet(dataset), dataset, cache}
{
}

PbiFilterQuery::PbiFilterQuery(const PbiFilter& filter, const DataSet& dataset)
    : PbiFilterQuery{filter, dataset, MakePbiIndexCache(dataset)}
{
}

PbiFilterQuery::PbiFilterQuery(const PbiFilter& filter, const DataSet& dataset,
                               const PbiIndexCache& cache)
    : internal::IQuery(), d_{std::make_unique<PbiFilterQueryPrivate>(filter, dataset, cache)}
{
}

PbiFilterQuery::~PbiFilterQuery() = default;

bool PbiFilterQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

uint32_t PbiFilterQuery::NumReads() const { return d_->reader_.NumReads(); }

}  // namespace BAM
}  // namespace PacBio
