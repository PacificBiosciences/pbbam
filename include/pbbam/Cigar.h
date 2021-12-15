#ifndef PBBAM_CIGAR_H
#define PBBAM_CIGAR_H

#include <pbbam/Config.h>

#include <pbbam/CigarOperation.h>

#include <pbcopper/data/Cigar.h>

namespace PacBio {
namespace BAM {

using Cigar PBBAM_DEPRECATED = PacBio::Data::Cigar;

PBBAM_DEPRECATED constexpr auto ReferenceLength = PacBio::Data::ReferenceLength;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_CIGAR_H
