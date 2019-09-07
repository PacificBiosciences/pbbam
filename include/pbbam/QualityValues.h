// File Description
/// \file QualityValues.h
/// \brief Defines the QualityValues class.
//
// Author: Derek Barnett

#ifndef QUALITYVALUES_H
#define QUALITYVALUES_H

#include "pbbam/Config.h"

#include <cstdint>

#include <string>
#include <vector>

#include <pbcopper/data/QualityValues.h>

#include "pbbam/QualityValue.h"

namespace PacBio {
namespace BAM {

using QualityValues PBBAM_DEPRECATED = PacBio::Data::QualityValues;

}  // namespace BAM
}  // namespace PacBio

#endif  // QUALITYVALUES_H
