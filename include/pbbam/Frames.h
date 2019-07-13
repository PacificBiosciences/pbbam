// File Description
/// \file Frames.h
/// \brief Defines the Frames class.
//
// Author: Derek Barnett

#ifndef FRAMES_H
#define FRAMES_H

#include "pbbam/Config.h"

#include <cstddef>
#include <cstdint>
#include <vector>

#include <pbcopper/data/Frames.h>

#ifndef PBBAM_NODEPRECATED_API

namespace PacBio {
namespace BAM {

using Frames = PacBio::Data::Frames;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_NODEPRECATED_API

#endif  // FRAMES_H
