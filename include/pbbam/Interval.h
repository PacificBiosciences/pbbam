#ifndef PBBAM_INTERVAL_H
#define PBBAM_INTERVAL_H

#include <pbbam/Config.h>

#include <pbcopper/data/Interval.h>

namespace PacBio {
namespace BAM {

using Interval PBBAM_DEPRECATED = Data::Interval;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_INTERVAL_H
