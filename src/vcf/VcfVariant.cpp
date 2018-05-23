// Author: Derek Barnett

#include <pbbam/vcf/VcfVariant.h>

#include <pbbam/vcf/VcfFormat.h>

namespace PacBio {
namespace VCF {

VcfVariant::VcfVariant(const std::string& text) { *this = VcfFormat::ParsedVariant(text); }

}  // namespace VCF
}  // namespace PacBio
