// File Description
/// \file CompositeFastaReader.inl
/// \brief Inline implementation for the composite FASTA reader, for
///        working with multiple input files.
//
// Author: Derek Barnett

#include "pbbam/CompositeFastaReader.h"

namespace PacBio {
namespace BAM {


inline CompositeFastaReader::CompositeFastaReader(const std::vector<std::string>& fastaFiles)
{
    for (auto&& fn : fastaFiles)
        readers_.emplace_back(new FastaReader{ fn });
}

inline CompositeFastaReader::CompositeFastaReader(std::vector<std::string>&& fastaFiles)
{
    for (auto&& fn : fastaFiles)
        readers_.emplace_back(new FastaReader{ std::move(fn) });
}

inline CompositeFastaReader::CompositeFastaReader(const DataSet& dataset)
    : CompositeFastaReader(dataset.FastaFiles())
{ }

inline bool CompositeFastaReader::GetNext(FastaSequence& seq)
{
    // try first reader, if successful return true
    // else pop reader and try next, until all readers exhausted
    while (!readers_.empty()) {
        auto& reader = readers_.front();
        if (reader->GetNext(seq))
            return true;
        else
            readers_.pop_front();
    }

    // no readers available
    return false;
}

} // namespace BAM
} // namespace PacBio
