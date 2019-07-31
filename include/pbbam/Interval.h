// File Description
/// \file Interval.h
/// \brief Defines the Interval class.
//
// Author: Derek Barnett

#ifndef INTERVAL_H
#define INTERVAL_H

#include "pbbam/Config.h"

#include <pbcopper/data/Interval.h>

namespace PacBio {
namespace BAM {

using Interval PBBAM_DEPRECATED = Data::Interval;

}  // namespace BAM
}  // namespace PacBio

#endif  // GENOMICINTERVAL_H
