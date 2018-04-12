// File Description
/// \file ClipType.h
/// \brief Defines the ClipType enum.
//
// Author: Derek Barnett

#ifndef CLIPTYPE_H
#define CLIPTYPE_H

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

#endif  // CLIPTYPE_H
