#ifndef PBBAM_CLIPTYPE_H
#define PBBAM_CLIPTYPE_H

#include <pbbam/Config.h>

namespace PacBio {
namespace BAM {

/// \brief This enum defines the modes supported by BamRecord clipping
///        operations.
///
/// Methods like BamRecord::Clip accept Position parameters - which may be in
/// either polymerase or reference coorindates. Using this enum as a flag
/// indicates how the positions should be interpreted.
///
enum class ClipType
{
    CLIP_NONE,         ///< No clipping will be performed.
    CLIP_TO_QUERY,     ///< Clipping positions are in polymerase coordinates.
    CLIP_TO_REFERENCE  ///< Clipping positions are in genomic coordinates.
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_CLIPTYPE_H
