// File Description
/// \file IFastqWriter.h
/// \brief Defines the IFastqWriter interface.
//
// Author: Derek Barnett

#ifndef IFASTQWRITER_H
#define IFASTQWRITER_H

#include "pbbam/Config.h"

#include <string>

#include <pbcopper/data/QualityValues.h>

#include "pbbam/IRecordWriter.h"

namespace PacBio {
namespace BAM {

class FastqSequence;

class IFastqWriter : public IRecordWriter
{
public:
    virtual ~IFastqWriter();

public:
    using IRecordWriter::Write;

    virtual void Write(const FastqSequence& fastq) = 0;
    virtual void Write(const std::string& name, const std::string& bases,
                       const Data::QualityValues& quals) = 0;
    virtual void Write(const std::string& name, const std::string& bases,
                       const std::string& quals) = 0;

protected:
    IFastqWriter();
};

}  // namespace BAM
}  // namespace PacBio

#endif  // IFASTQWRITER_H
