// File Description
/// \file QualityValue.h
/// \brief Defines the QualityValue class.
//
// Author: Derek Barnett

#ifndef QUALITYVALUE_H
#define QUALITYVALUE_H

#include "pbbam/Config.h"

#include <cstdint>

#include <string>
#include <vector>

#include <pbcopper/data/QualityValue.h>

namespace PacBio {
namespace BAM {

using QualityValue PBBAM_DEPRECATED = PacBio::Data::QualityValue;

}  // namespace BAM
}  // namespace PacBio

#endif  // QUALITYVALUE_H
