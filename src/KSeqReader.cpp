// File Description
/// \file KSeqReader.cpp
/// \brief Implements the KSeqReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "KSeqReader.h"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace PacBio {
namespace BAM {

KSeqReader::KSeqReader(const std::string& fn)
    : fp_{gzopen(fn.c_str(), "r")}, seq_{kseq_init(fp_.get())}
{
    // check file handle
    if (fp_ == nullptr) {
        std::ostringstream msg;
        msg << "[pbbam] kseq FASTX reader ERROR: could not open file:\n"
            << "  file: " << fn << '\n';
        throw std::runtime_error{msg.str()};
    }

    // check kseq sequence handle
    assert(seq_ != nullptr);
}

KSeqReader::KSeqReader(KSeqReader&&) noexcept = default;

KSeqReader& KSeqReader::operator=(KSeqReader&&) noexcept = default;

KSeqReader::~KSeqReader() = default;

std::string KSeqReader::Bases() const { return std::string{seq_->seq.s, seq_->seq.l}; }

std::string KSeqReader::Name() const { return std::string{seq_->name.s, seq_->name.l}; }

std::string KSeqReader::Qualities() const { return std::string{seq_->qual.s, seq_->qual.l}; }

bool KSeqReader::ReadNext()
{
    const auto result = kseq_read(seq_.get());
    if (result == -1)  // EOF
        return false;
    return true;
}

}  // namespace BAM
}  // namespace PacBio
