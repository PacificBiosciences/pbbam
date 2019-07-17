// File Description
/// \file CigarOperation.h
/// \brief Defines the CigarOperationType enum & CigarOperation class.
//
// Author: Derek Barnett

#ifndef CIGAROPERATION_H
#define CIGAROPERATION_H

#include "pbbam/Config.h"

#include <cstdint>
#include <stdexcept>

#include <pbcopper/data/CigarOperation.h>

namespace PacBio {
namespace BAM {

using CigarOperation PBBAM_DEPRECATED = PacBio::Data::CigarOperation;
using CigarOperationType PBBAM_DEPRECATED = PacBio::Data::CigarOperationType;

PBBAM_DEPRECATED constexpr auto ConsumesQuery = PacBio::Data::ConsumesQuery;
PBBAM_DEPRECATED constexpr auto ConsumesReference = PacBio::Data::ConsumesReference;

}  // namespace BAM
}  // namespace PacBio

#endif  // CIGAROPERATION_H
