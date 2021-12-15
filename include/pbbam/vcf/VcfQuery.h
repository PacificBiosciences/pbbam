#ifndef PBBAM_VCF_VCFQUERY_H
#define PBBAM_VCF_VCFQUERY_H

#include <pbbam/Config.h>

#include <pbbam/internal/QueryBase.h>
#include <pbbam/vcf/VcfFile.h>
#include <pbbam/vcf/VcfReader.h>
#include <pbbam/vcf/VcfVariant.h>

#include <string>

namespace PacBio {
namespace VCF {

class VcfQuery : public BAM::internal::QueryBase<VcfVariant>
{
public:
    explicit VcfQuery(std::string fn);
    explicit VcfQuery(const VcfFile& file);

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(VcfVariant& var) override;

private:
    VcfReader reader_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFQUERY_H
