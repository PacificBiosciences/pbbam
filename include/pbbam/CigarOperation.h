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

using CigarOperation = PacBio::Data::CigarOperation;
using CigarOperationType = PacBio::Data::CigarOperationType;

constexpr auto ConsumesQuery = PacBio::Data::ConsumesQuery;
constexpr auto ConsumesReference = PacBio::Data::ConsumesReference;

}  // namespace BAM
}  // namespace PacBio

#endif  // CIGAROPERATION_H
