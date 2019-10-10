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

namespace PacBio {
namespace BAM {

using Frames PBBAM_DEPRECATED = PacBio::Data::Frames;

}  // namespace BAM
}  // namespace PacBio

#endif  // FRAMES_H
