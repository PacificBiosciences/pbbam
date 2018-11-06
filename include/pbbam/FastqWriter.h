// File Description
/// \file FastqWriter.h
/// \brief Defines the FastqWriter class.
//
// Author: Derek Barnett

#ifndef FASTQWRITER_H
#define FASTQWRITER_H

#include <fstream>
#include <iostream>
#include <string>

namespace PacBio {
namespace BAM {

class BamRecord;
class FastqSequence;
class QualityValues;

class FastqWriter
{
public:
    FastqWriter(const std::string& fn);

public:
    void Write(const FastqSequence& fastq);

    void Write(const BamRecord& bam);

    void Write(const std::string& name, const std::string& bases, const QualityValues& quals);

    void Write(const std::string& name, const std::string& bases, const std::string& quals);

private:
    std::ofstream file_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTQWRITER_H
