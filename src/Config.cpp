// File Description
/// \file Config.cpp
/// \brief Initializes global variable defaults.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

// Initialized to -1 to indicate default. Client code may set this or not.
//
// To respect client code or else fallback to default[OFF], this value should be used like this:
//
//    hts_verbose = ( PacBio::BAM::HtslibVerbosity == -1 ? 0 : PacBio::BAM::HtslibVerbosity);
//
//
//
int HtslibVerbosity = -1;

}  // namespace BAM
}  // namespace PacBio
