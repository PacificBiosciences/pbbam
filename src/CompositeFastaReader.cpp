// File Description
/// \file BamRecordView.cpp
/// \brief Implements the BamRecordTags utility class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/CompositeFastaReader.h"

namespace PacBio {
namespace BAM {

CompositeFastaReader::CompositeFastaReader(const std::vector<std::string>& fastaFiles)
{
    for (const auto& fn : fastaFiles)
        readers_.emplace_back(std::make_unique<FastaReader>(fn));
}

CompositeFastaReader::CompositeFastaReader(const DataSet& dataset)
    : CompositeFastaReader{dataset.FastaFiles()}
{
}

bool CompositeFastaReader::GetNext(FastaSequence& seq)
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

}  // namespace BAM
}  // namespace PacBio
