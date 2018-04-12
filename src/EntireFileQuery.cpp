// File Description
/// \file EntireFileQuery.cpp
/// \brief Implements the EntireFileQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/EntireFileQuery.h"

#include "pbbam/CompositeBamReader.h"

namespace PacBio {
namespace BAM {

struct EntireFileQuery::EntireFileQueryPrivate
{
    EntireFileQueryPrivate(const DataSet &dataset) : reader_(dataset) {}

    SequentialCompositeBamReader reader_;
};

EntireFileQuery::EntireFileQuery(const DataSet &dataset)
    : internal::IQuery(), d_(new EntireFileQueryPrivate(dataset))
{
}

EntireFileQuery::~EntireFileQuery() {}

bool EntireFileQuery::GetNext(BamRecord &r) { return d_->reader_.GetNext(r); }

}  // namespace BAM
}  // namespace PacBio
