#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfQuery.h>

#include <type_traits>

namespace PacBio {
namespace VCF {

VcfQuery::VcfQuery(std::string fn) : VcfQuery{VcfFile{std::move(fn)}} {}

VcfQuery::VcfQuery(const VcfFile& file) : BAM::internal::QueryBase<VcfVariant>(), reader_{file} {}

bool VcfQuery::GetNext(VcfVariant& var) { return reader_.GetNext(var); }

}  // namespace VCF
}  // namespace PacBio
