#include "PbbamInternalConfig.h"

#include <pbbam/EntireFileQuery.h>

#include <pbbam/CompositeBamReader.h>

namespace PacBio {
namespace BAM {

class EntireFileQuery::EntireFileQueryPrivate
{
public:
    EntireFileQueryPrivate(const DataSet& dataset) : reader_{dataset} {}

    SequentialCompositeBamReader reader_;
};

EntireFileQuery::EntireFileQuery(const DataSet& dataset)
    : internal::IQuery{}, d_(new EntireFileQueryPrivate(dataset))
{}

EntireFileQuery::EntireFileQuery(EntireFileQuery&&) noexcept = default;

EntireFileQuery& EntireFileQuery::operator=(EntireFileQuery&&) noexcept = default;

EntireFileQuery::~EntireFileQuery() = default;

bool EntireFileQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

}  // namespace BAM
}  // namespace PacBio
