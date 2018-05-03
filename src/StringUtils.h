// Author: Derek Barnett

#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include "pbbam/SamTagCodec.h"
#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {
namespace internal {

inline std::string MakeSamTag(std::string tag, std::string value)
{
    return PacBio::BAM::MakeSamTag(std::move(tag), std::move(value));
}

inline std::vector<std::string> Split(const std::string& line, const char delim = '\t')
{
    return PacBio::BAM::Split(line, delim);
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // STRINGUTILS_H
