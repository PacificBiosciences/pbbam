// File Description
/// \file GenomicInterval.h
/// \brief Defines the GenomicInterval class.
//
// Author: Derek Barnett

#ifndef GENOMICINTERVAL_H
#define GENOMICINTERVAL_H

#include "pbbam/Config.h"

#include <pbcopper/data/GenomicInterval.h>

namespace PacBio {
namespace BAM {

using GenomicInterval PBBAM_DEPRECATED = Data::GenomicInterval;

}  // namespace BAM
}  // namespace PacBio

#endif  // GENOMICINTERVAL_H
