// File Description
/// \file BarcodeQuery.cpp
/// \brief Implements the BarcodeQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BarcodeQuery.h"

#include <cstdint>

#include "pbbam/CompositeBamReader.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiFilterTypes.h"

namespace PacBio {
namespace BAM {

struct BarcodeQuery::BarcodeQueryPrivate
{
    BarcodeQueryPrivate(const int16_t barcode, const DataSet& dataset)
        : reader_{PbiBarcodeFilter{barcode}, dataset}
    {
    }

    PbiFilterCompositeBamReader<Compare::None> reader_;  // unsorted
};

BarcodeQuery::BarcodeQuery(const int16_t barcode, const DataSet& dataset)
    : internal::IQuery(), d_{std::make_unique<BarcodeQueryPrivate>(barcode, dataset)}
{
}

BarcodeQuery::~BarcodeQuery() {}

bool BarcodeQuery::GetNext(BamRecord& r) { return d_->reader_.GetNext(r); }

}  // namespace BAM
}  // namespace PacBio
