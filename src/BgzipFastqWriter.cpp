#include "PbbamInternalConfig.h"

#include <pbbam/BgzipFastqWriter.h>

#include <pbbam/BamRecord.h>
#include <pbbam/FastqSequence.h>
#include <pbbam/FormatUtils.h>

#include <sstream>
#include <stdexcept>

namespace PacBio {
namespace BAM {

BgzipFastqWriter::BgzipFastqWriter(const std::string& fn)
    : BgzipFastqWriter{fn, BgzipWriterConfig{}}
{
}

BgzipFastqWriter::BgzipFastqWriter(const std::string& fn, const BgzipWriterConfig& config)
    : IFastqWriter{}, writer_{fn, config}
{
    if (!FormatUtils::IsFastqFilename(fn)) {
        std::ostringstream s;
        s << "[pbbam] bgzipped FASTQ writer ERROR: not a recognized FASTQ extension:\n"
          << "  file: " << fn;
        throw std::runtime_error{s.str()};
    }
}

void BgzipFastqWriter::Write(const FastqSequence& fastq)
{
    Write(fastq.Name(), fastq.Bases(), fastq.Qualities());
}

void BgzipFastqWriter::Write(const BamRecord& bam)
{
    Write(bam.FullName(), bam.Sequence(), bam.Qualities());
}

void BgzipFastqWriter::Write(const BamRecordImpl& bam)
{
    Write(bam.Name(), bam.Sequence(), bam.Qualities());
}

void BgzipFastqWriter::Write(const std::string& name, const std::string& bases,
                             const Data::QualityValues& quals)
{
    Write(name, bases, quals.Fastq());
}

void BgzipFastqWriter::Write(const std::string& name, const std::string& bases,
                             const std::string& quals)
{
    std::string out{"@" + name + '\n' + bases + "\n+\n"};
    if (!quals.empty()) {
        out += quals;
    } else {
        out += std::string(bases.size(), '!');
    }
    out.push_back('\n');
    writer_.Write(out);
}

}  // namespace BAM
}  // namespace PacBio
