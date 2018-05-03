// File Description
/// \file FastaSequenceQuery.cpp
/// \brief Implements the FastaSequenceQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastaSequenceQuery.h"

#include "pbbam/CompositeFastaReader.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

struct FastaSequenceQuery::FastaSequenceQueryPrivate
{
    FastaSequenceQueryPrivate(const DataSet& dataset) : reader_{dataset} {}

    CompositeFastaReader reader_;
};

FastaSequenceQuery::FastaSequenceQuery(const DataSet& dataset)
    : internal::QueryBase<FastaSequence>(), d_{std::make_unique<FastaSequenceQueryPrivate>(dataset)}
{
}

FastaSequenceQuery::~FastaSequenceQuery() {}

bool FastaSequenceQuery::GetNext(FastaSequence& seq) { return d_->reader_.GetNext(seq); }

}  // namespace BAM
}  // namespace PacBio
