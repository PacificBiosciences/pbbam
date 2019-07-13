// File Description
/// \file Strand.h
/// \brief Defines the Strand enum.
//
// Author: Derek Barnett

#ifndef STRAND_H
#define STRAND_H

#include "pbbam/Config.h"

#include <pbcopper/data/Strand.h>

#ifndef PBBAM_NODEPRECATED_API

namespace PacBio {
namespace BAM {

using Strand = PacBio::Data::Strand;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_NODEPRECATED_API
#endif  // STRAND_H
