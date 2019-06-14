// File Description
/// \file FastqWriter.h
/// \brief Defines the FastqWriter class.
//
// Author: Derek Barnett

#ifndef FASTQWRITER_H
#define FASTQWRITER_H

#include "pbbam/Config.h"

#include <fstream>
#include <iostream>
#include <string>

#include "pbbam/IFastqWriter.h"

namespace PacBio {
namespace BAM {

class BamRecord;
class BamRecordImpl;
class FastqSequence;
class QualityValues;

class FastqWriter final : public IFastqWriter
{
public:
    FastqWriter(const std::string& fn);

public:
    // IFastqWriter
    void Write(const FastqSequence& fastq);
    void Write(const std::string& name, const std::string& bases, const QualityValues& quals);
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

#endif  // FASTQWRITER_H
