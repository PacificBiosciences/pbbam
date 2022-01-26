#ifndef PBBAM_IFASTAWRITER_H
#define PBBAM_IFASTAWRITER_H

#include <pbbam/Config.h>

#include <pbbam/FastaSequence.h>
#include <pbbam/IRecordWriter.h>

#include <string>

namespace PacBio {
namespace BAM {

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

#endif  // PBBAM_IFASTAWRITER_H
