// File Description
/// \file GenomicIntervalQuery.cpp
/// \brief Implements the GenomicIntervalQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/GenomicIntervalQuery.h"

#include "pbbam/CompositeBamReader.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

struct GenomicIntervalQuery::GenomicIntervalQueryPrivate
{
    GenomicIntervalQueryPrivate(const GenomicInterval& interval, const DataSet& dataset)
        : reader_{interval, dataset}
    {
    }

    GenomicIntervalCompositeBamReader reader_;
};

GenomicIntervalQuery::GenomicIntervalQuery(const GenomicInterval& interval, const DataSet& dataset)
    : internal::IQuery(), d_{std::make_unique<GenomicIntervalQueryPrivate>(interval, dataset)}
{
}

GenomicIntervalQuery::~GenomicIntervalQuery() {}

bool GenomicIntervalQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

GenomicIntervalQuery& GenomicIntervalQuery::Interval(const GenomicInterval& interval)
{
    d_->reader_.Interval(interval);
    return *this;
}

const GenomicInterval& GenomicIntervalQuery::Interval() const { return d_->reader_.Interval(); }

}  // namespace BAM
}  // namespace PacBio
