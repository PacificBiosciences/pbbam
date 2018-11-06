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

namespace PacBio {
namespace BAM {

class BamRecord;
class FastaSequence;

class FastaWriter
{
public:
    FastaWriter(const std::string& fn);

public:
    void Write(const FastaSequence& fastq);

    void Write(const BamRecord& bam);

    void Write(const std::string& name, const std::string& bases);

private:
    std::ofstream file_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTAWRITER_H
