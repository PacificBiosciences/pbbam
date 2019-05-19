// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFQUERY_H
#define PBBAM_VCF_VCFQUERY_H

#include "pbbam/Config.h"

#include <string>

#include <pbbam/vcf/VcfFile.h>
#include <pbbam/vcf/VcfReader.h>
#include <pbbam/vcf/VcfVariant.h>

#include <pbbam/internal/QueryBase.h>

namespace PacBio {
namespace VCF {

class VcfQuery : public PacBio::BAM::internal::QueryBase<VcfVariant>
{
public:
    explicit VcfQuery(std::string fn);
    explicit VcfQuery(const VcfFile& file);

    VcfQuery() = delete;
    VcfQuery(const VcfQuery&) = delete;
    VcfQuery(VcfQuery&&) noexcept;
    VcfQuery& operator=(const VcfQuery&) = delete;
    VcfQuery& operator=(VcfQuery&&) noexcept(std::is_nothrow_move_assignable<VcfReader>::value);
    ~VcfQuery();

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
