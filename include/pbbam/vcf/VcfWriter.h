// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFWRITER_H
#define PBBAM_VCF_VCFWRITER_H

#include <memory>
#include <string>

namespace PacBio {
namespace VCF {

class VcfHeader;
class VcfVariant;

class VcfWriter
{
public:
    VcfWriter(std::string filename, const VcfHeader& header);

    VcfWriter() = delete;
    VcfWriter(const VcfWriter&) = delete;
    VcfWriter(VcfWriter&&) = default;
    VcfWriter& operator=(const VcfWriter&) = delete;
    VcfWriter& operator=(VcfWriter&&) = default;
    ~VcfWriter();

public:
    bool Write(const VcfVariant& var);

private:
    struct VcfWriterPrivate;
    std::unique_ptr<VcfWriterPrivate> d_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFWRITER_H
