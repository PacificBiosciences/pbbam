#include "PbbamInternalConfig.h"

#include <pbbam/Config.h>

#include <stdexcept>
#include <string>
#include <tuple>

#include <htslib/hts.h>
#include <pbcopper/data/CigarOperation.h>

#include <pbbam/StringUtilities.h>

namespace PacBio {
namespace BAM {

// Disable htslib's own logging at startup. Client code can still override with
// hts_set_log_level(HTS_LOG_FOO).
static const int DisableHtslibLogging = []() {
    hts_set_log_level(HTS_LOG_OFF);
    return 0;
}();

bool DoesHtslibSupportLongCigar()
{
    const std::string htsVersion = hts_version();

    // remove any "-<blah>" for non-release versions
    const auto versionBase = PacBio::BAM::Split(htsVersion, '-');
    if (versionBase.empty()) {
        throw std::runtime_error{"[pbbam] config ERROR: invalid htslib version format: '" +
                                 htsVersion + "'"};
    }

    // grab major/minor version numbers
    const auto versionParts = PacBio::BAM::Split(versionBase[0], '.');
    if (versionParts.size() < 2) {
        throw std::runtime_error{"[pbbam] config ERROR: invalid htslib version format: '" +
                                 htsVersion + "'"};
    }

    // check against v1.7
    const int versionMajor = std::stoi(versionParts[0]);
    const int versionMinor = std::stoi(versionParts[1]);
    static constexpr const int v17_major = 1;
    static constexpr const int v17_minor = 7;
    return std::tie(versionMajor, versionMinor) >= std::tie(v17_major, v17_minor);
}

#ifdef PBBAM_PERMISSIVE_CIGAR
static const bool PermissiveCigar = []() {
    Data::CigarOperation::DisableAutoValidation();
    return true;
}();
#endif  // PBBAM_PERMISSIVE_CIGAR

}  // namespace BAM
}  // namespace PacBio
