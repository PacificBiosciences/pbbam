// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFFILE_H
#define PBBAM_VCF_VCFFILE_H

#include "pbbam/Config.h"

#include <string>

#include <pbbam/vcf/VcfHeader.h>

namespace PacBio {
namespace VCF {

class VcfFile
{
public:
    explicit VcfFile(std::string fn);

public:
    const std::string& Filename() const;
    const VcfHeader& Header() const;

private:
    std::string filename_;
    VcfHeader header_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFFILE_H
