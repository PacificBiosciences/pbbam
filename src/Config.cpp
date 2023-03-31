#include "PbbamInternalConfig.h"

#include <pbbam/Config.h>

#include <pbbam/StringUtilities.h>

#include <pbcopper/data/CigarOperation.h>

#include <htslib/hts.h>

#include <stdexcept>
#include <string>
#include <tuple>

namespace PacBio {
namespace BAM {

// Disable htslib's own logging at startup. Client code can still override with
// hts_set_log_level(HTS_LOG_FOO).
static const int DisableHtslibLogging = []() {
    hts_set_log_level(HTS_LOG_OFF);
    return 0;
}();

#ifdef PBBAM_PERMISSIVE_CIGAR
static const bool PermissiveCigar = []() {
    Data::CigarOperation::DisableAutoValidation();
    return true;
}();
#endif  // PBBAM_PERMISSIVE_CIGAR

}  // namespace BAM
}  // namespace PacBio
