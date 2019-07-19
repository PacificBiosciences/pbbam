// File Description
/// \file Accuracy.h
/// \brief Defines the Accuracy class.
//
// Author: Derek Barnett

#ifndef ACCURACY_H
#define ACCURACY_H

#include "pbbam/Config.h"

#include <pbcopper/data/Accuracy.h>

namespace PacBio {
namespace BAM {

using Accuracy PBBAM_DEPRECATED = Data::Accuracy;

}  // namespace BAM
}  // namespace PacBio

#endif  // ACCURACY_H
