// File Description
/// \file Cigar.h
/// \brief Defines the Cigar class.
//
// Author: Derek Barnett

#ifndef CIGAR_H
#define CIGAR_H

#include "pbbam/Config.h"

#include <string>
#include <vector>

#include <pbcopper/data/Cigar.h>
#include <pbcopper/data/CigarOperation.h>

#ifndef PBBAM_NODEPRECATED_API

namespace PacBio {
namespace BAM {

using Cigar = PacBio::Data::Cigar;

constexpr auto ReferenceLength = PacBio::Data::ReferenceLength;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_NODEPRECATED_API

#endif  // CIGAR_H
