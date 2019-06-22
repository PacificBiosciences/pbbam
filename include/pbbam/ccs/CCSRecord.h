// Author: Derek Barnett

#ifndef PBBAM_CCS_CCSRECORD_H
#define PBBAM_CCS_CCSRECORD_H

#include "pbbam/Config.h"

#include "pbbam/Accuracy.h"
#include "pbbam/Frames.h"
#include "pbbam/LocalContextFlags.h"
#include "pbbam/Position.h"
#include "pbbam/SNR.h"

namespace PacBio {
namespace CCS {

struct CCSRecord
{
    int32_t HoleNumber = 0;

    PacBio::BAM::Position QueryStart = 0;

    PacBio::BAM::Position QueryEnd = 0;

    PacBio::BAM::LocalContextFlags LocalContextFlags =
        PacBio::BAM::LocalContextFlags::NO_LOCAL_CONTEXT;

    PacBio::BAM::Accuracy Accuracy = 0.0f;

    PacBio::BAM::SNR SignalToNoise = {0.0, 0.0, 0.0, 0.0};

    std::string Sequence;

    PacBio::BAM::Frames PulseWidths;
};

}  // namespace CCS
}  // namespace PacBio

#endif  // PBBAM_CCS_CCSRECORD_H
