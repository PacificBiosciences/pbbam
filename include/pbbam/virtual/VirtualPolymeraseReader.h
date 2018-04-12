// File Description
/// \file VirtualPolymeraseReader.h
/// \brief Defines the VirtualPolymeraseReader class.
//
// Author: Armin Töpfer

#ifndef VIRTUALPOLYMERASEREADER_H
#define VIRTUALPOLYMERASEREADER_H

#include "pbbam/virtual/VirtualPolymeraseBamRecord.h"
#include "pbbam/virtual/ZmwReadStitcher.h"

namespace PacBio {
namespace BAM {

/// \deprecated Use ZmwReadStitcher instead.
typedef ZmwReadStitcher VirtualPolymeraseReader;

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALPOLYMERASEREADER_H
