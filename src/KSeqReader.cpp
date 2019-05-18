// File Description
/// \file KSeqReader.cpp
/// \brief Implements the KSeqReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "KSeqReader.h"

#include <stdexcept>

namespace PacBio {
namespace BAM {

KSeqReader::KSeqReader(const std::string& fn) : fp_{gzopen(fn.c_str(), "r")}, seq_{kseq_init(fp_)}
{
    if (fp_ == nullptr || seq_.get() == nullptr)
        throw std::runtime_error{"KSeqReader: could not open file for reading: " + fn};
}

KSeqReader::KSeqReader(KSeqReader&&) = default;

KSeqReader& KSeqReader::operator=(KSeqReader&&) = default;

KSeqReader::~KSeqReader() { gzclose(fp_); }

bool KSeqReader::ReadNext()
{
    const auto result = kseq_read(seq_.get());
    if (result == -1)  // EOF
        return false;
    return true;
}

}  // namespace BAM
}  // namespace PacBio
