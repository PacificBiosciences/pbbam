#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfQuery.h>

#include <cassert>

#include <type_traits>

namespace PacBio {
namespace VCF {

static_assert(!std::is_copy_constructible<VcfQuery>::value,
              "VcfQuery(const VcfQuery&) is not = delete");
static_assert(!std::is_copy_assignable<VcfQuery>::value,
              "VcfQuery& operator=(const VcfQuery&) is not = delete");

static_assert(std::is_nothrow_move_constructible<VcfQuery>::value ==
                  std::is_nothrow_move_constructible<VcfReader>::value,
              "");
static_assert(std::is_nothrow_move_assignable<VcfQuery>::value ==
                  std::is_nothrow_move_assignable<VcfReader>::value,
              "");

VcfQuery::VcfQuery(std::string fn) : VcfQuery{VcfFile{std::move(fn)}} {}

VcfQuery::VcfQuery(const VcfFile& file)
    : PacBio::BAM::internal::QueryBase<VcfVariant>(), reader_{file}
{
}

bool VcfQuery::GetNext(VcfVariant& var) { return reader_.GetNext(var); }

}  // namespace VCF
}  // namespace PacBio
