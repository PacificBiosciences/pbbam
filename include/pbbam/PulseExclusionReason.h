// File Description
/// \file PulseExclusionReason.h
/// \brief Defines the PulseExclusionReason enum.
//
// Author: Derek Barnett

#ifndef PULSE_EXCLUSION_REASON_H
#define PULSE_EXCLUSION_REASON_H

#include "pbbam/Config.h"

#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief This enum defines the possible pulse exclusion reasons
///
enum class PulseExclusionReason : uint8_t
{
    BASE = 0,
    SHORT_PULSE,
    BURST,
    PAUSE
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PULSE_EXCLUSION_REASON_H
