// File Description
/// \file BgzipFastaWriter.cpp
/// \brief Implements the BgzipFastaWriter class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BgzipFastaWriter.h"

#include <stdexcept>

#include "pbbam/BamRecord.h"
#include "pbbam/FastqSequence.h"
#include "pbbam/FormatUtils.h"

namespace PacBio {
namespace BAM {

BgzipFastaWriter::BgzipFastaWriter(const std::string& fn) : IFastaWriter{}, writer_{fn}
{
    if (!FormatUtils::IsFastaFilename(fn)) {
        throw std::runtime_error{"BgzipFastaWriter: filename '" + fn +
                                 "' is not recognized as a FASTA file."};
    }
}

BgzipFastaWriter::BgzipFastaWriter(const std::string& fn, const BgzipWriterConfig& config)
    : IFastaWriter{}, writer_{fn, config}
{
    if (!FormatUtils::IsFastaFilename(fn)) {
        throw std::runtime_error{"BgzipFastaWriter: filename '" + fn +
                                 "' is not recognized as a FASTA file."};
    }
}

void BgzipFastaWriter::TryFlush() {}

void BgzipFastaWriter::Write(const BamRecordImpl& bam) { Write(bam.Name(), bam.Sequence()); }

void BgzipFastaWriter::Write(const FastaSequence& fastq) { Write(fastq.Name(), fastq.Bases()); }

void BgzipFastaWriter::Write(const BamRecord& bam) { Write(bam.FullName(), bam.Sequence()); }

void BgzipFastaWriter::Write(const std::string& name, const std::string& bases)
{
    // TODO: wrap bases
    std::string out{">" + name + '\n' + bases + '\n'};
    writer_.Write(out);
}

}  // namespace BAM
}  // namespace PacBio
