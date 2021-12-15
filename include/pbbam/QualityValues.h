#ifndef PBBAM_QUALITYVALUES_H
#define PBBAM_QUALITYVALUES_H

#include <pbbam/Config.h>

#include <pbbam/QualityValue.h>

#include <pbcopper/data/QualityValues.h>

namespace PacBio {
namespace BAM {

using QualityValues PBBAM_DEPRECATED = PacBio::Data::QualityValues;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_QUALITYVALUES_H
