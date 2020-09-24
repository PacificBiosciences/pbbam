#include "PbbamInternalConfig.h"

#include <pbbam/vcf/VcfFile.h>

#include <cassert>

#include <type_traits>

#include <pbbam/vcf/VcfFormat.h>

namespace PacBio {
namespace VCF {

static_assert(std::is_copy_constructible<VcfFile>::value,
              "VcfFile(const VcfFile&) is not = default");
static_assert(std::is_copy_assignable<VcfFile>::value,
              "VcfFile& operator=(const VcfFile&) is not = default");

static_assert(std::is_nothrow_move_constructible<VcfFile>::value,
              "VcfFile(VcfFile&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<VcfFile>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

VcfFile::VcfFile(std::string fn)
    : filename_{std::move(fn)}, header_{VcfFormat::HeaderFromFile(filename_)}
{
}

const std::string& VcfFile::Filename() const { return filename_; }

const VcfHeader& VcfFile::Header() const { return header_; }

}  // namespace VCF
}  // namespace PacBio
