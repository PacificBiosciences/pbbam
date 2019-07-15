// Author: Derek Barnett

#ifndef PBBAM_CCS_CCSRECORD_H
#define PBBAM_CCS_CCSRECORD_H

#include "pbbam/Config.h"

#include <pbcopper/data/Frames.h>
#include <pbcopper/data/Position.h>
#include <pbcopper/data/SNR.h>

#include "pbbam/Accuracy.h"
#include "pbbam/LocalContextFlags.h"

namespace PacBio {
namespace CCS {

struct CCSRecord
{
    int32_t HoleNumber = 0;

    Data::Position QueryStart = 0;

    Data::Position QueryEnd = 0;

    PacBio::BAM::LocalContextFlags LocalContextFlags =
        PacBio::BAM::LocalContextFlags::NO_LOCAL_CONTEXT;

    PacBio::BAM::Accuracy Accuracy = 0.0f;

    Data::SNR SignalToNoise = {0.0, 0.0, 0.0, 0.0};

    std::string Sequence;

    Data::Frames PulseWidths;
};

}  // namespace CCS
}  // namespace PacBio

#endif  // PBBAM_CCS_CCSRECORD_H
