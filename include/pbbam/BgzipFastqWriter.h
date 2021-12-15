#ifndef PBBAM_BGZIPFASTQWRITER_H
#define PBBAM_BGZIPFASTQWRITER_H

#include <pbbam/Config.h>

#include <pbbam/BgzipWriter.h>
#include <pbbam/IFastqWriter.h>

#include <string>

namespace PacBio {
namespace BAM {

class BgzipFastqWriter final : public IFastqWriter
{
public:
    explicit BgzipFastqWriter(const std::string& fn);
    BgzipFastqWriter(const std::string& fn, const BgzipWriterConfig& config);

public:
    // IFastqWriter
    void Write(const FastqSequence& fastq);
    void Write(const std::string& name, const std::string& bases, const Data::QualityValues& quals);
    void Write(const std::string& name, const std::string& bases, const std::string& quals);

    // IRecordWriter
    void Write(const BamRecord& bam);
    void Write(const BamRecordImpl& bam);

private:
    BgzipWriter writer_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BGZFFASTQWRITER_H
