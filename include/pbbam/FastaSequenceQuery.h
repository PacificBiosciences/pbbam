// File Description
/// \file FastaSequenceQuery.h
/// \brief Defines the FastaSequenceQuery class.
//
// Author: Derek Barnett

#ifndef FASTASEQUENCEQUERY_H
#define FASTASEQUENCEQUERY_H

#include "pbbam/DataSet.h"
#include "pbbam/FastaSequence.h"
#include "pbbam/internal/QueryBase.h"

#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

///
/// \brief The FastaSequence class represents a FASTA record (name & bases)
///
class FastaSequenceQuery : public internal::QueryBase<FastaSequence>
{
public:
    /// \name Constructors & Related Methods
    /// \{

    FastaSequenceQuery(const PacBio::BAM::DataSet& dataset);
    ~FastaSequenceQuery() override;

    /// \}

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(FastaSequence& seq) override;

private:
    struct FastaSequenceQueryPrivate;
    std::unique_ptr<FastaSequenceQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTASEQUENCEQUERY_H
