#ifndef PBBAM_CIGAROPERATION_H
#define PBBAM_CIGAROPERATION_H

#include <pbbam/Config.h>

#include <pbcopper/data/CigarOperation.h>

namespace PacBio {
namespace BAM {

using CigarOperation PBBAM_DEPRECATED = PacBio::Data::CigarOperation;
using CigarOperationType PBBAM_DEPRECATED = PacBio::Data::CigarOperationType;

PBBAM_DEPRECATED constexpr auto ConsumesQuery = PacBio::Data::ConsumesQuery;
PBBAM_DEPRECATED constexpr auto ConsumesReference = PacBio::Data::ConsumesReference;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_CIGAROPERATION_H
