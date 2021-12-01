#include "PbbamInternalConfig.h"

#include <pbbam/GenomicIntervalQuery.h>

#include <pbbam/CompositeBamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/GenomicInterval.h>

namespace PacBio {
namespace BAM {

class GenomicIntervalQuery::GenomicIntervalQueryPrivate
{
public:
    GenomicIntervalQueryPrivate(const DataSet& dataset, const BaiIndexCache& cache)
        : reader_{dataset, cache}
    {
    }

    GenomicIntervalQueryPrivate(const Data::GenomicInterval& interval, const DataSet& dataset,
                                const BaiIndexCache& cache)
        : reader_{interval, dataset, cache}
    {
    }

    GenomicIntervalCompositeBamReader reader_;
};

GenomicIntervalQuery::GenomicIntervalQuery(const DataSet& dataset)
    : GenomicIntervalQuery{dataset, MakeBaiIndexCache(dataset)}
{
}

GenomicIntervalQuery::GenomicIntervalQuery(const DataSet& dataset, const BaiIndexCache& cache)
    : internal::IQuery{}, d_{std::make_unique<GenomicIntervalQueryPrivate>(dataset, cache)}
{
}

GenomicIntervalQuery::GenomicIntervalQuery(const Data::GenomicInterval& interval,
                                           const DataSet& dataset)
    : GenomicIntervalQuery{interval, dataset, MakeBaiIndexCache(dataset)}
{
}

GenomicIntervalQuery::GenomicIntervalQuery(const Data::GenomicInterval& interval,
                                           const DataSet& dataset, const BaiIndexCache& cache)
    : internal::IQuery{}
    , d_{std::make_unique<GenomicIntervalQueryPrivate>(interval, dataset, cache)}
{
}

GenomicIntervalQuery::GenomicIntervalQuery(GenomicIntervalQuery&&) noexcept = default;

GenomicIntervalQuery& GenomicIntervalQuery::operator=(GenomicIntervalQuery&&) noexcept = default;

GenomicIntervalQuery::~GenomicIntervalQuery() = default;

bool GenomicIntervalQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

GenomicIntervalQuery& GenomicIntervalQuery::Interval(const Data::GenomicInterval& interval)
{
    d_->reader_.Interval(interval);
    return *this;
}

const Data::GenomicInterval& GenomicIntervalQuery::Interval() const
{
    return d_->reader_.Interval();
}

}  // namespace BAM
}  // namespace PacBio
