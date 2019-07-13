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

#ifndef PBBAM_NODEPRECATED_API

namespace PacBio {
namespace BAM {

using QualityValue = PacBio::Data::QualityValue;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_NODEPRECATED_API

#endif  // QUALITYVALUE_H
