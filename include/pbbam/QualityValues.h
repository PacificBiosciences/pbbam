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

#ifndef PBBAM_NODEPRECATED_API

namespace PacBio {
namespace BAM {

using QualityValues = PacBio::Data::QualityValues;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_NODEPRECATED_API

#endif  // QUALITYVALUES_H
