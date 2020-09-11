#ifndef PBBAM_STRAND_H
#define PBBAM_STRAND_H

#include <pbbam/Config.h>

#include <pbcopper/data/Strand.h>

namespace PacBio {
namespace BAM {

using Strand PBBAM_DEPRECATED = PacBio::Data::Strand;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_STRAND_H
