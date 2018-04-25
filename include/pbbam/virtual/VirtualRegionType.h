// File Description
/// \file VirtualRegionType.h
/// \brief Defines the VirtualRegionType enum.
//
// Author: Derek Barnett

#ifndef REGIONTYPE_H
#define REGIONTYPE_H

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief This enum defines the types of annotated region.
///
enum class VirtualRegionType  // : char
{
    ADAPTER = 0x41,   ///< Adapter region ('A')
    BARCODE = 0x42,   ///< Barcode region ('B')
    FILTERED = 0x46,  ///< Filtered subread ('F')
    SUBREAD = 0x53,   ///< Subread ('S')
    HQREGION = 0x48,  ///< High-quality region ('H')
    LQREGION = 0x4C   ///< Low-quality region ('L'), i.e. outside the HQ region
};

}  // namespace BAM
}  // namespace PacBio

#endif  // REGIONTYPE_H
