#ifndef PBBAM_VIRTUALPOLYMERASECOMPOSITEREADER_H
#define PBBAM_VIRTUALPOLYMERASECOMPOSITEREADER_H

#include <pbbam/Config.h>

#include <pbbam/virtual/ZmwReadStitcher.h>

namespace PacBio {
namespace BAM {

/// \deprecated Use ZmwReadStitcher instead.
using VirtualPolymeraseCompositeReader = ZmwReadStitcher;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_VIRTUALPOLYMERASECOMPOSITEREADER_H
