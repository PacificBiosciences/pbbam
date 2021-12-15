#include "PbbamInternalConfig.h"

#include <pbbam/vcf/VcfFile.h>

#include <pbbam/vcf/VcfFormat.h>

#include <type_traits>

namespace PacBio {
namespace VCF {

VcfFile::VcfFile(std::string fn)
    : filename_{std::move(fn)}, header_{VcfFormat::HeaderFromFile(filename_)}
{
}

const std::string& VcfFile::Filename() const { return filename_; }

const VcfHeader& VcfFile::Header() const { return header_; }

}  // namespace VCF
}  // namespace PacBio
