// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFSORT_H
#define PBBAM_VCF_VCFSORT_H

#include <string>

#include "pbbam/Config.h"
#include "pbbam/vcf/VcfFile.h"

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
inline void SortFile(const std::string& inputFilename, const std::string& outputFilename)
{
    SortFile(VcfFile{inputFilename}, outputFilename);
}

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFSORT_H
