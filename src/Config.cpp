// File Description
/// \file Config.cpp
/// \brief Initializes global variable defaults.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Config.h"

#include <stdexcept>
#include <string>
#include <tuple>

#include <htslib/hts.h>

#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {

// Initialized to -1 to indicate default. Client code may set this or not.
//
// To respect client code or else fallback to default[OFF], this value should be used like this:
//
//    hts_verbose = ( PacBio::BAM::HtslibVerbosity == -1 ? 0 : PacBio::BAM::HtslibVerbosity);
//
//
//
int HtslibVerbosity = -1;

bool DoesHtslibSupportLongCigar()
{
    const std::string htsVersion = hts_version();

    // remove any "-<blah>" for non-release versions
    const auto versionBase = PacBio::BAM::Split(htsVersion, '-');
    if (versionBase.empty())
        throw std::runtime_error{"invalid htslib version format: " + htsVersion};

    // grab major/minor version numbers
    const auto versionParts = PacBio::BAM::Split(versionBase[0], '.');
    if (versionParts.size() < 2)
        throw std::runtime_error{"invalid htslib version format: " + htsVersion};

    // check against v1.7
    const int versionMajor = std::stoi(versionParts[0]);
    const int versionMinor = std::stoi(versionParts[1]);
    static constexpr const int v17_major = 1;
    static constexpr const int v17_minor = 7;
    return std::tie(versionMajor, versionMinor) >= std::tie(v17_major, v17_minor);
}

}  // namespace BAM
}  // namespace PacBio
