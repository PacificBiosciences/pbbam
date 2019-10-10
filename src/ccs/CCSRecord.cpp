// File Description
/// \file CCSRecord.cpp
/// \brief Implements the CCSRecordclass.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ccs/CCSRecord.h"

namespace PacBio {
namespace CCS {

// clang-format off
Data::Read CCSRecord::ToRead(std::string movieName, std::string chemistry) const
{
    return Data::Read{
        Data::ReadId{movieName, static_cast<size_t>(HoleNumber),
                     Data::Interval{QueryStart, QueryEnd}},
        Sequence,
        PulseWidths,
        boost::none,
        LocalContextFlags,
        Accuracy,
        SignalToNoise,
        std::move(chemistry)};
}
// clang-format on

}  // namespace CCS
}  // namespace PacBio
