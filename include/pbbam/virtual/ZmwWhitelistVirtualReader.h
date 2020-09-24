#ifndef PBBAM_ZMWWHITELISTVIRTUALREADER_H
#define PBBAM_ZMWWHITELISTVIRTUALREADER_H

#include <pbbam/Config.h>

#include <pbbam/virtual/WhitelistedZmwReadStitcher.h>

namespace PacBio {
namespace BAM {

/// \deprecated Use WhitelistedZmwReadStitcher instead.
using ZmwWhitelistVirtualReader = WhitelistedZmwReadStitcher;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWWHITELISTVIRTUALREADER_H
