// Author: Derek Barnett

#include <pbbam/vcf/VcfHeader.h>

#include <pbbam/vcf/VcfFormat.h>

namespace PacBio {
namespace VCF {

VcfHeader::VcfHeader() { Version(VcfFormat::CurrentVersion()); }

VcfHeader::VcfHeader(const std::string& hdrText) { *this = VcfFormat::ParsedHeader(hdrText); }

}  // namespace VCF
}  // namespace PacBio
