// File Description
/// \file VirtualPolymeraseCompositeReader.h
/// \brief Defines the VirtualPolymeraseCompositeReader class.
//
// Author: Derek Barnett

#ifndef VIRTUALPOLYMERASECOMPOSITEREADER_H
#define VIRTUALPOLYMERASECOMPOSITEREADER_H

#include "pbbam/virtual/VirtualPolymeraseBamRecord.h"
#include "pbbam/virtual/ZmwReadStitcher.h"

namespace PacBio {
namespace BAM {

/// \deprecated Use ZmwReadStitcher instead.
typedef ZmwReadStitcher VirtualPolymeraseCompositeReader;

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALPOLYMERASECOMPOSITEREADER_H
