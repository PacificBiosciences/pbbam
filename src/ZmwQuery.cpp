#include "PbbamInternalConfig.h"

#include <pbbam/ZmwQuery.h>

#include <cstdint>

#include <pbbam/CompositeBamReader.h>
#include <pbbam/PbiFilterTypes.h>

namespace PacBio {
namespace BAM {

class ZmwQuery::ZmwQueryPrivate
{
public:
    ZmwQueryPrivate(std::vector<int32_t> zmwWhitelist, const DataSet& dataset)
        : reader_{PbiZmwFilter{std::move(zmwWhitelist)}, dataset}
    {
    }

    PbiFilterCompositeBamReader<Compare::Zmw> reader_;
};

ZmwQuery::ZmwQuery(std::vector<int32_t> zmwWhitelist, const DataSet& dataset)
    : internal::IQuery(), d_{std::make_unique<ZmwQueryPrivate>(zmwWhitelist, dataset)}
{
}

ZmwQuery::~ZmwQuery() = default;

bool ZmwQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

}  // namespace BAM
}  // namespace PacBio
