// File Description
/// \file FastqWriter.cpp
/// \brief Implements the FastqWriter class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqWriter.h"

#include "pbbam/BamRecord.h"
#include "pbbam/FastqSequence.h"
#include "pbbam/QualityValues.h"

namespace PacBio {
namespace BAM {

FastqWriter::FastqWriter(const std::string& fn) : IRecordWriter(), file_{fn}
{
    if (!file_) throw std::runtime_error{"FastqWriter could not open for writing: " + fn};
}

void FastqWriter::TryFlush() { file_.flush(); }

void FastqWriter::Write(const FastqSequence& fastq)
{
    Write(fastq.Name(), fastq.Bases(), fastq.Qualities());
}

void FastqWriter::Write(const BamRecord& bam)
{
    Write(bam.FullName(), bam.Sequence(), bam.Qualities());
}

void FastqWriter::Write(const BamRecordImpl& bam)
{
    Write(bam.Name(), bam.Sequence(), bam.Qualities());
}

void FastqWriter::Write(const std::string& name, const std::string& bases,
                        const QualityValues& quals)
{
    Write(name, bases, quals.Fastq());
}

void FastqWriter::Write(const std::string& name, const std::string& bases, const std::string& quals)
{

    file_ << "@" << name << '\n' << bases << '\n' << "+\n";

    if (!quals.empty())
        file_ << quals;
    else {
        std::string q(bases.size(), '!');
        file_ << q;
    }

    file_ << '\n';
}

}  // namespace BAM
}  // namespace PacBio
