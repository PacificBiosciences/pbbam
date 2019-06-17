// File Description
/// \file IFastaWriter.h
/// \brief Defines the IFastaWriter interface.
//
// Author: Derek Barnett

#ifndef IFASTAWRITER_H
#define IFASTAWRITER_H

#include "pbbam/Config.h"

#include <string>

#include "pbbam/IRecordWriter.h"

namespace PacBio {
namespace BAM {

class FastaSequence;

class IFastaWriter : public IRecordWriter
{
public:
    virtual ~IFastaWriter();

public:
    using IRecordWriter::Write;

    virtual void Write(const FastaSequence& fastq) = 0;
    virtual void Write(const std::string& name, const std::string& bases) = 0;

protected:
    IFastaWriter();
};

}  // namespace BAM
}  // namespace PacBio

#endif  // IFASTAWRITER_H
