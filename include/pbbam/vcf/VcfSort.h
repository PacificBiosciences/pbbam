#ifndef PBBAM_VCF_VCFSORT_H
#define PBBAM_VCF_VCFSORT_H

#include <pbbam/Config.h>

#include <pbbam/vcf/VcfFile.h>

#include <string>

namespace PacBio {
namespace VCF {

///
/// \brief SortFile
/// \param file
/// \param outputFilename
///
void SortFile(const VcfFile& file, const std::string& outputFilename);

///
/// \brief SortFile
/// \param inputFilename
/// \param outputFilename
///
void SortFile(const std::string& inputFilename, const std::string& outputFilename);

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFSORT_H
