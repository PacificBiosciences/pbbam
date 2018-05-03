// File Description
/// \file SubreadLengthQuery.cpp
/// \brief Implements the SubreadLengthQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SubreadLengthQuery.h"

#include <cstdint>

#include "pbbam/CompositeBamReader.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiFilterTypes.h"

namespace PacBio {
namespace BAM {

struct SubreadLengthQuery::SubreadLengthQueryPrivate
{
    SubreadLengthQueryPrivate(const int32_t length, const Compare::Type compareType,
                              const DataSet& dataset)
        : reader_(PbiQueryLengthFilter(length, compareType), dataset)
    {
    }

    PbiFilterCompositeBamReader<Compare::None> reader_;  // unsorted
};

SubreadLengthQuery::SubreadLengthQuery(const int32_t length, const Compare::Type compareType,
                                       const DataSet& dataset)
    : internal::IQuery()
    , d_{std::make_unique<SubreadLengthQueryPrivate>(length, compareType, dataset)}
{
}

SubreadLengthQuery::~SubreadLengthQuery() {}

bool SubreadLengthQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

uint32_t SubreadLengthQuery::NumReads() const { return d_->reader_.NumReads(); }

}  // namespace BAM
}  // namespace PacBio
