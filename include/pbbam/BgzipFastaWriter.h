#ifndef PBBAM_BGZIPFASTAWRITER_H
#define PBBAM_BGZIPFASTAWRITER_H

#include <pbbam/Config.h>

#include <pbbam/BgzipWriter.h>
#include <pbbam/IFastaWriter.h>

#include <string>

namespace PacBio {
namespace BAM {

class BgzipFastaWriter final : public IFastaWriter
{
public:
    explicit BgzipFastaWriter(const std::string& fn);
    BgzipFastaWriter(const std::string& fn, const BgzipWriterConfig& config);

public:
    // IFastaWriter
    void Write(const FastaSequence& fastq);
    void Write(const std::string& name, const std::string& bases);

    // IRecordWriter
    void Write(const BamRecord& bam);
    void Write(const BamRecordImpl& bam);

private:
    BgzipWriter writer_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BGZIPFASTAWRITER_H
