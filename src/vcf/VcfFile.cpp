#include "PbbamInternalConfig.h"

#include "pbbam/vcf/VcfFile.h"

#include "pbbam/vcf/VcfFormat.h"

namespace PacBio {
namespace VCF {

VcfFile::VcfFile(std::string fn)
    : filename_{std::move(fn)}, header_{VcfFormat::HeaderFromFile(filename_)}
{
}

VcfFile::VcfFile(const VcfFile&) = default;

VcfFile::VcfFile(VcfFile&&) noexcept = default;

VcfFile& VcfFile::operator=(const VcfFile&) = default;

VcfFile& VcfFile::operator=(VcfFile&&) PBBAM_NOEXCEPT_MOVE_ASSIGN = default;

VcfFile::~VcfFile() = default;

const std::string& VcfFile::Filename() const { return filename_; }

const VcfHeader& VcfFile::Header() const { return header_; }

}  // namespace VCF
}  // namespace PacBio
