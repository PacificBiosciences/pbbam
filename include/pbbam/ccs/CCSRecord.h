#ifndef PBBAM_CCS_CCSRECORD_H
#define PBBAM_CCS_CCSRECORD_H

#include <pbbam/Config.h>

#include <pbcopper/data/Accuracy.h>
#include <pbcopper/data/Frames.h>
#include <pbcopper/data/LocalContextFlags.h>
#include <pbcopper/data/Position.h>
#include <pbcopper/data/Read.h>
#include <pbcopper/data/SNR.h>

namespace PacBio {
namespace CCS {

struct CCSRecord
{
    int32_t HoleNumber = 0;

    Data::Position QueryStart = 0;

    Data::Position QueryEnd = 0;

    Data::LocalContextFlags LocalContextFlags = Data::LocalContextFlags::NO_LOCAL_CONTEXT;

    Data::Accuracy Accuracy = 0.0f;

    Data::SNR SignalToNoise = {0.0, 0.0, 0.0, 0.0};

    std::string Sequence;

    Data::Frames PulseWidths;

    ///
    /// Create a Data::Read from this CCSRecord.
    ///
    Data::Read ToRead(std::string movieName, std::string chemistry) const;
};

}  // namespace CCS
}  // namespace PacBio

#endif  // PBBAM_CCS_CCSRECORD_H
