// File Description
/// \file VirtualPolymeraseBamRecord.h
/// \brief Defines the VirtualPolymeraseBamRecord class.
//
// Author: Armin Töpfer

#ifndef VIRTUALPOLYMERASEBAMRECORD_H
#define VIRTUALPOLYMERASEBAMRECORD_H

#include "pbbam/virtual/VirtualZmwBamRecord.h"

namespace PacBio {
namespace BAM {

/// \deprecated Use VirtualZmwBamRecord instead.
using VirtualPolymeraseBamRecord = VirtualZmwBamRecord;

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALPOLYMERASEBAMRECORD_H
