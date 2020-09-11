#ifndef PBBAM_GENOMICINTERVAL_H
#define PBBAM_GENOMICINTERVAL_H

#include <pbbam/Config.h>

#include <pbcopper/data/GenomicInterval.h>

namespace PacBio {
namespace BAM {

using GenomicInterval PBBAM_DEPRECATED = Data::GenomicInterval;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_GENOMICINTERVAL_H
