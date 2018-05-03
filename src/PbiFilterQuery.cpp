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

struct PbiFilterQuery::PbiFilterQueryPrivate
{
    PbiFilterQueryPrivate(const PbiFilter& filter, const DataSet& dataset)
        : reader_{filter, dataset}
    {
    }

    PbiFilterCompositeBamReader<Compare::None> reader_;  // unsorted
};

PbiFilterQuery::PbiFilterQuery(const DataSet& dataset)
    : PbiFilterQuery{PbiFilter::FromDataSet(dataset), dataset}
{
}

PbiFilterQuery::PbiFilterQuery(const PbiFilter& filter, const DataSet& dataset)
    : internal::IQuery(), d_{std::make_unique<PbiFilterQueryPrivate>(filter, dataset)}
{
}

PbiFilterQuery::~PbiFilterQuery() {}

bool PbiFilterQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

uint32_t PbiFilterQuery::NumReads() const { return d_->reader_.NumReads(); }

}  // namespace BAM
}  // namespace PacBio
