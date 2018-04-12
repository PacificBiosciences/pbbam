// File Description
/// \file PulseBehavior.h
/// \brief Defines the PulseBehavior enum.
//
// Author: Derek Barnett

#ifndef PULSEBEHAVIOR_H
#define PULSEBEHAVIOR_H

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

#endif  // PULSEBEHAVIOR_H
