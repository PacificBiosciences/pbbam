// File Description
/// \file BgzipFastqWriter.h
/// \brief Defines the BgzipFastqWriter class.
//
// Author: Derek Barnett

#ifndef BGZIPFASTQWRITER_H
#define BGZIPFASTQWRITER_H

#include "pbbam/Config.h"

#include <fstream>
#include <iostream>
#include <string>

#include "pbbam/BgzipWriter.h"
#include "pbbam/IRecordWriter.h"

namespace PacBio {
namespace BAM {

class BamRecord;
class BamRecordImpl;
class FastqSequence;
class QualityValues;

class BgzipFastqWriter : public IRecordWriter
{
public:
    explicit BgzipFastqWriter(const std::string& fn);
    BgzipFastqWriter(const std::string& fn, const BgzipWriterConfig& config);

public:
    void Write(const FastqSequence& fastq);
    void Write(const std::string& name, const std::string& bases, const QualityValues& quals);
    void Write(const std::string& name, const std::string& bases, const std::string& quals);

    // IRecordWriter
    void TryFlush() override;
    void Write(const BamRecord& bam) override;
    void Write(const BamRecordImpl& bam) override;

private:
    BgzipWriter writer_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BGZFFASTQWRITER_H
