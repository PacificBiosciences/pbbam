// File Description
/// \file FastaWriter.cpp
/// \brief Implements the FastaWriter class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastaWriter.h"

#include <sstream>
#include <stdexcept>

#include "pbbam/BamRecord.h"
#include "pbbam/FastqSequence.h"
#include "pbbam/FormatUtils.h"

namespace PacBio {
namespace BAM {

FastaWriter::FastaWriter(const std::string& fn) : IFastaWriter{}
{
    if (!FormatUtils::IsFastaFilename(fn)) {
        std::ostringstream s;
        s << "[pbbam] FASTA writer ERROR: not a recognized FASTA extension:\n"
          << "  file: " << fn;
        throw std::runtime_error{s.str()};
    }

    file_.open(fn);
    if (!file_) {
        std::ostringstream s;
        s << "[pbbam] FASTA writer ERROR: could not open file for writing:\n"
          << "  file: " << fn;
        throw std::runtime_error{s.str()};
    }
}

void FastaWriter::TryFlush() { file_.flush(); }

void FastaWriter::Write(const BamRecordImpl& bam) { Write(bam.Name(), bam.Sequence()); }

void FastaWriter::Write(const FastaSequence& fastq) { Write(fastq.Name(), fastq.Bases()); }

void FastaWriter::Write(const BamRecord& bam) { Write(bam.FullName(), bam.Sequence()); }

void FastaWriter::Write(const std::string& name, const std::string& bases)
{
    // TODO: wrap bases
    file_ << ">" << name << '\n' << bases << '\n';
}

}  // namespace BAM
}  // namespace PacBio
