#include "PbbamInternalConfig.h"

#include <pbbam/FastqWriter.h>

#include <pbbam/BamRecord.h>
#include <pbbam/FastqSequence.h>
#include <pbbam/FormatUtils.h>
#include <pbbam/QualityValues.h>
#include "ErrnoReason.h"

#include <sstream>
#include <stdexcept>

namespace PacBio {
namespace BAM {

FastqWriter::FastqWriter(const std::string& fn) : IFastqWriter{}
{
    if (!FormatUtils::IsFastqFilename(fn)) {
        std::ostringstream s;
        s << "[pbbam] FASTQ writer ERROR: not a recognized FASTQ extension:\n"
          << "  file: " << fn;
        throw std::runtime_error{s.str()};
    }

    file_.open(fn);
    if (!file_) {
        std::ostringstream s;
        s << "[pbbam] FASTQ writer ERROR: could not open file for writing:\n"
          << "  file: " << fn;
        MaybePrintErrnoReason(s);
        throw std::runtime_error{s.str()};
    }
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
                        const Data::QualityValues& quals)
{
    Write(name, bases, quals.Fastq());
}

void FastqWriter::Write(const std::string& name, const std::string& bases, const std::string& quals)
{
    file_ << "@" << name << '\n' << bases << '\n' << "+\n";

    if (!quals.empty()) {
        file_ << quals;
    } else {
        std::string q(bases.size(), '!');
        file_ << q;
    }

    file_ << '\n';
}

}  // namespace BAM
}  // namespace PacBio
