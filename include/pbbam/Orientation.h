// File Description
/// \file Orientation.h
/// \brief Defines the Orientation enum.
//
// Author: Derek Barnett

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include "pbbam/Config.h"

#include <pbcopper/data/Orientation.h>

namespace PacBio {
namespace BAM {

using Orientation PBBAM_DEPRECATED = PacBio::Data::Orientation;

}  // namespace BAM
}  // namespace PacBio

#endif  // ORIENTATION_H
