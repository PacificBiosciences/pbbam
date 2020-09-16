#ifndef PBBAM_VIRTUALPOLYMERASEREADER_H
#define PBBAM_VIRTUALPOLYMERASEREADER_H

#include <pbbam/Config.h>

#include <pbbam/virtual/ZmwReadStitcher.h>

namespace PacBio {
namespace BAM {

/// \deprecated Use ZmwReadStitcher instead.
using VirtualPolymeraseReader = ZmwReadStitcher;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_VIRTUALPOLYMERASEREADER_H
