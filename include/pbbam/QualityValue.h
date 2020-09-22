#ifndef PBBAM_QUALITYVALUE_H
#define PBBAM_QUALITYVALUE_H

#include <pbbam/Config.h>

#include <pbcopper/data/QualityValue.h>

namespace PacBio {
namespace BAM {

using QualityValue PBBAM_DEPRECATED = PacBio::Data::QualityValue;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_QUALITYVALUE_H
