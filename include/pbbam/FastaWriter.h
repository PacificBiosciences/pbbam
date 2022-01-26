#ifndef PBBAM_FASTAWRITER_H
#define PBBAM_FASTAWRITER_H

#include <pbbam/Config.h>

#include <pbbam/IFastaWriter.h>

#include <fstream>
#include <string>

namespace PacBio {
namespace BAM {

class FastaWriter final : public IFastaWriter
{
public:
    FastaWriter(const std::string& fn);

public:
    // IFastaWriter
    void Write(const FastaSequence& fastq);
    void Write(const std::string& name, const std::string& bases);

    // IRecordWriter
    void TryFlush();
    void Write(const BamRecord& bam);
    void Write(const BamRecordImpl& bam);

private:
    std::ofstream file_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FASTAWRITER_H
