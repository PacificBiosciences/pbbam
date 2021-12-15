#ifndef PBBAM_FASTQWRITER_H
#define PBBAM_FASTQWRITER_H

#include <pbbam/Config.h>

#include <pbbam/IFastqWriter.h>

#include <fstream>
#include <string>

namespace PacBio {
namespace BAM {

class FastqWriter final : public IFastqWriter
{
public:
    FastqWriter(const std::string& fn);

public:
    // IFastqWriter
    void Write(const FastqSequence& fastq);
    void Write(const std::string& name, const std::string& bases, const Data::QualityValues& quals);
    void Write(const std::string& name, const std::string& bases, const std::string& quals);

    // IRecordWriter
    void TryFlush();
    void Write(const BamRecord& bam);
    void Write(const BamRecordImpl& bam);

private:
    std::ofstream file_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FASTQWRITER_H
