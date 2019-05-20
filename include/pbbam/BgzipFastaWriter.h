// File Description
/// \file BgzipFastaWriter.h
/// \brief Defines the BgzipFastaWriter class.
//
// Author: Derek Barnett

#ifndef BGZIPFASTAWRITER_H
#define BGZIPFASTAWRITER_H

#include "pbbam/Config.h"

#include <string>

#include "pbbam/BgzipWriter.h"
#include "pbbam/IRecordWriter.h"

namespace PacBio {
namespace BAM {

class BamRecord;
class BamRecordImpl;
class FastaSequence;

class BgzipFastaWriter : public IRecordWriter
{
public:
    explicit BgzipFastaWriter(const std::string& fn);
    BgzipFastaWriter(const std::string& fn, const BgzipWriterConfig& config);

public:
    void Write(const FastaSequence& fastq);
    void Write(const std::string& name, const std::string& bases);

    // IRecordWriter
    void TryFlush() override;
    void Write(const BamRecord& bam) override;
    void Write(const BamRecordImpl& bam) override;

private:
    BgzipWriter writer_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BGZIPFASTAWRITER_H
