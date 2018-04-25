// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "FofnReader.h"

#include <iostream>

namespace PacBio {
namespace BAM {
namespace internal {

std::vector<std::string> FofnReader::Files(std::istream& in)
{
    std::vector<std::string> files;
    std::string fn;
    while (std::getline(in, fn))
        files.push_back(fn);
    return files;
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
