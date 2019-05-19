// File Description
/// \file VirtualPolymeraseReader.h
/// \brief Defines the VirtualPolymeraseReader class.
//
// Author: Armin Töpfer

#ifndef VIRTUALPOLYMERASEREADER_H
#define VIRTUALPOLYMERASEREADER_H

#include "pbbam/Config.h"

#include "pbbam/virtual/ZmwReadStitcher.h"

namespace PacBio {
namespace BAM {

/// \deprecated Use ZmwReadStitcher instead.
using VirtualPolymeraseReader = ZmwReadStitcher;

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALPOLYMERASEREADER_H
