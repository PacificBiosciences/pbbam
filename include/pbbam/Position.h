// File Description
/// \file Position.h
/// \brief Defines the Position typedef.
//
// Author: Derek Barnett

#ifndef POSITION_H
#define POSITION_H

#include <cstdint>
#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief This type is used to refer to genomic positions.
/// \typedef typedef int32_t PacBio::BAM::Position
///
/// We use a signed integer because SAM/BAM uses the -1 value to indicate
/// unknown or unmapped positions.
///
using Position = int32_t;

/// \brief This constant is widely used as a "missing" or "invalid" position
///        marker.
///
static const Position UnmappedPosition{-1};

}  // namespace BAM
}  // namespace PacBio

#endif  // POSITION_H
