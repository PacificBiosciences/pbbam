#ifndef PBBAM_VCF_VCFFILE_INL
#define PBBAM_VCF_VCFFILE_INL

#include <pbbam/vcf/VcfFile.h>
#include <pbbam/vcf/VcfFormat.h>

namespace PacBio {
namespace VCF {

inline VcfFile::VcfFile(std::string fn)
    : filename_{std::move(fn)}
    , header_{VcfFormat::HeaderFromFile(filename_)}
{ }

inline const std::string& VcfFile::Filename() const { return filename_; }

inline const VcfHeader& VcfFile::Header() const { return header_; }

} // namespace VCF
} // namespace PacBio

#endif // PBBAM_VCF_VCFFILE_INL
