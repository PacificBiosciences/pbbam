// File Description
/// \file FastaWriter.h
/// \brief Defines the FastaWriter class.
//
// Author: Derek Barnett

#ifndef FASTAWRITER_H
#define FASTAWRITER_H

#include <fstream>
#include <iostream>
#include <string>

#include "pbbam/IRecordWriter.h"

namespace PacBio {
namespace BAM {

class BamRecord;
class BamRecordImpl;
class FastaSequence;

class FastaWriter : public IRecordWriter
{
public:
    FastaWriter(const std::string& fn);

public:
    void Write(const FastaSequence& fastq);
    void Write(const std::string& name, const std::string& bases);

    // IRecordWriter
    void TryFlush() override;
    void Write(const BamRecord& bam) override;
    void Write(const BamRecordImpl& bam) override;

private:
    std::ofstream file_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTAWRITER_H
