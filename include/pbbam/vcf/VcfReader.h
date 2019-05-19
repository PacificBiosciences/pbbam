// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFREADER_H
#define PBBAM_VCF_VCFREADER_H

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "pbbam/Config.h"

#include <pbbam/vcf/VcfFile.h>
#include <pbbam/vcf/VcfHeader.h>
#include <pbbam/vcf/VcfVariant.h>

namespace PacBio {
namespace VCF {

///
/// \brief The VcfReader class
///
class VcfReader
{
public:
    explicit VcfReader(std::string fn);
    explicit VcfReader(const VcfFile& file);

    VcfReader() = delete;
    VcfReader(const VcfReader&) = delete;
    VcfReader(VcfReader&&) noexcept;
    VcfReader& operator=(const VcfReader&) = delete;
    VcfReader& operator=(VcfReader&&) PBBAM_NOEXCEPT_MOVE_ASSIGN;
    ~VcfReader();

public:
    const VcfHeader& Header() const;

    bool GetNext(VcfVariant& var);

private:
    void FetchNext();

private:
    std::unique_ptr<std::ifstream> in_;
    VcfHeader header_;
    std::string line_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFREADER_H
