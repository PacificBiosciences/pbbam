#ifndef PBBAM_VCF_VCFREADER_H
#define PBBAM_VCF_VCFREADER_H

#include <pbbam/Config.h>

#include <fstream>
#include <memory>
#include <string>

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

public:
    const VcfHeader& Header() const;

    bool GetNext(VcfVariant& var);

private:
    void FetchNext();

private:
    std::ifstream in_;
    VcfHeader header_;
    std::string line_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFREADER_H
