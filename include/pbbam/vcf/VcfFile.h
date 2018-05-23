// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFFILE_H
#define PBBAM_VCF_VCFFILE_H

#include <string>

#include <pbbam/vcf/VcfHeader.h>

namespace PacBio {
namespace VCF {

class VcfFile
{
public:
    explicit VcfFile(std::string fn);

    VcfFile() = delete;
    VcfFile(const VcfFile&) = default;
    VcfFile(VcfFile&&) = default;
    VcfFile& operator=(const VcfFile&) = default;
    VcfFile& operator=(VcfFile&&) = default;
    ~VcfFile() = default;

public:
    const std::string& Filename() const;
    const VcfHeader& Header() const;

private:
    std::string filename_;
    VcfHeader header_;
};

}  // namespace VCF
}  // namespace PacBio

#include "pbbam/vcf/internal/VcfFile.inl"

#endif  // PBBAM_VCF_VCFFILE_H
