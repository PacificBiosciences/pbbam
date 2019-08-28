// File Description
/// \file FastaWriter.h
/// \brief Defines the FastaWriter class.
//
// Author: Derek Barnett

#ifndef FASTAWRITER_H
#define FASTAWRITER_H

#include "pbbam/Config.h"

#include <fstream>
#include <iostream>
#include <string>

#include "pbbam/IFastaWriter.h"

namespace PacBio {
namespace BAM {

class BamRecord;
class BamRecordImpl;
class FastaSequence;

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

#endif  // FASTAWRITER_H
