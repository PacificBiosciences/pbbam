// File Description
/// \file CompositeFastaReader.h
/// \brief Defines the composite FASTA reader, for working with multiple input
///       files.
//
// Author: Derek Barnett

#ifndef COMPOSITEFASTAREADER_H
#define COMPOSITEFASTAREADER_H

#include "pbbam/Config.h"
#include "pbbam/DataSet.h"
#include "pbbam/FastaReader.h"
#include "pbbam/FastaSequence.h"

#include <deque>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The CompositeFastaReader class provides read access to
///        multiple FASTA files, reading through the entire contents of each
///        file.
///
/// Input files will be accessed in the order provided to the constructor. Each
/// file's contents will be exhausted before moving on to the next one (as
/// opposed to a "round-robin" scheme).
///
class PBBAM_EXPORT CompositeFastaReader
{
public:
    /// \name Contstructors & Related Methods
    /// \{

    CompositeFastaReader(const std::vector<std::string>& fastaFiles);
    CompositeFastaReader(const DataSet& dataset);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Fetches next FASTA sequence.
    ///
    /// \returns true on success, false if no more data available.
    ///
    bool GetNext(FastaSequence& seq);

    /// \}

private:
    std::deque<std::unique_ptr<FastaReader> > readers_;
};

}  // namespace BAM
}  // namespace PacBio

#include "internal/CompositeFastaReader.inl"

#endif  // COMPOSITEFASTAREADER_H
