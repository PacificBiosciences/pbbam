#ifndef PBBAM_PBBAMVERSION_H
#define PBBAM_PBBAMVERSION_H

#include <pbbam/Config.h>

#include <string>
#include <tuple>

namespace PacBio {
namespace BAM {

std::string LibraryGitSha1String();
std::string LibraryVersionString();
std::tuple<int, int, int> LibraryVersionTriple();

inline std::string LibraryFormattedVersion()
{
    return LibraryVersionString() + " (commit " + LibraryGitSha1String() + ")";
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PBBAMVERSION_H
