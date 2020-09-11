#ifndef PBBAM_PULSEBEHAVIOR_H
#define PBBAM_PULSEBEHAVIOR_H

#include <pbbam/Config.h>

namespace PacBio {
namespace BAM {

/// \brief This enum defines the pulsecall modes supported by BamRecord tag
///        accessors.
///
enum class PulseBehavior
{
    BASECALLS_ONLY,  ///< "Squashed" pulses not included, only basecalls.
    ALL              ///< All pulses included.
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PULSEBEHAVIOR_H
