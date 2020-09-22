#ifndef PBBAM_SNR_H
#define PBBAM_SNR_H

#include <pbbam/Config.h>

#include <pbcopper/data/SNR.h>

namespace PacBio {
namespace BAM {

using SNR PBBAM_DEPRECATED = PacBio::Data::SNR;

PBBAM_DEPRECATED constexpr auto ClampSNR = PacBio::Data::ClampSNR;

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_SNR_H
