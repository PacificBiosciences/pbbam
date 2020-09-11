#ifndef PBBAM_IFASTQWRITER_H
#define PBBAM_IFASTQWRITER_H

#include <pbbam/Config.h>

#include <string>

#include <pbbam/FastqSequence.h>
#include <pbbam/IRecordWriter.h>

namespace PacBio {
namespace BAM {

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

#endif  // PBBAM_IFASTQWRITER_H
