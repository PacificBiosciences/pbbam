#ifndef PBBAM_ZMWTYPE_H
#define PBBAM_ZMWTYPE_H

#include <pbbam/Config.h>

namespace PacBio {
namespace BAM {

/// \brief This enum defines the different ZMW categories of scraps
///
enum class ZmwType : char
{
    CONTROL = 'C',
    MALFORMED = 'M',
    NORMAL = 'N',
    SENTINEL = 'S'
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWTYPE_H
