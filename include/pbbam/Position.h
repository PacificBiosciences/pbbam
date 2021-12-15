#ifndef PBBAM_POSITION_H
#define PBBAM_POSITION_H

#include <pbbam/Config.h>

#include <pbcopper/data/Position.h>

#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief This type is used to refer to genomic positions.
/// \typedef typedef int32_t PacBio::BAM::Position
///
/// We use a signed integer because SAM/BAM uses the -1 value to indicate
/// unknown or unmapped positions.
///
using Position PBBAM_DEPRECATED = PacBio::Data::Position;

/// \brief This constant is widely used as a "missing" or "invalid" position
///        marker.
///
PBBAM_DEPRECATED constexpr Position UnmappedPosition{-1};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_POSITION_H
