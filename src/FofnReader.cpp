#include "PbbamInternalConfig.h"

#include "FofnReader.h"

#include <istream>

namespace PacBio {
namespace BAM {

std::vector<std::string> FofnReader::Files(std::istream& in)
{
    std::vector<std::string> files;
    std::string fn;
    while (std::getline(in, fn))
        files.push_back(fn);
    return files;
}

}  // namespace BAM
}  // namespace PacBio
