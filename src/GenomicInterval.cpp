// File Description
/// \file GenomicInterval.cpp
/// \brief Implements the GenomicInterval class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/GenomicInterval.h"

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "StringUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

// returns sequence name & sets begin/end, from input regionString
std::string parseRegionString(const std::string& reg, PacBio::BAM::Position* begin,
                              PacBio::BAM::Position* end)
{
    const std::vector<std::string> parts = internal::Split(reg, ':');
    if (parts.empty() || parts.size() > 2) throw std::runtime_error{"malformed region string"};

    // given name only, default min,max intervals
    if (parts.size() == 1) {
        *begin = 0;
        *end = 1 << 29;
    }

    // parse interval from input
    else if (parts.size() == 2) {
        const std::vector<std::string> intervalParts = internal::Split(parts.at(1), '-');
        if (intervalParts.empty() || intervalParts.size() > 2)
            throw std::runtime_error{"malformed region string"};
        *begin = std::stoi(intervalParts.at(0));
        *end = std::stoi(intervalParts.at(1));
    }

    return parts.at(0);
}

}  // namespace internal

GenomicInterval::GenomicInterval(std::string name, Position start, Position stop)
    : name_{std::move(name)}, interval_{std::move(start), std::move(stop)}
{
}

GenomicInterval::GenomicInterval(const std::string& samtoolsRegionString)
{
    Position begin;
    Position end;
    name_ = internal::parseRegionString(samtoolsRegionString, &begin, &end);
    interval_ = PacBio::BAM::Interval<Position>(begin, end);
}

bool GenomicInterval::CoveredBy(const GenomicInterval& other) const
{
    if (name_ != other.name_) return false;
    return interval_.CoveredBy(other.interval_);
}

bool GenomicInterval::Covers(const GenomicInterval& other) const
{
    if (name_ != other.name_) return false;
    return interval_.Covers(other.interval_);
}

bool GenomicInterval::Intersects(const GenomicInterval& other) const
{
    if (name_ != other.name_) return false;
    return interval_.Intersects(other.interval_);
}

}  // namespace BAM
}  // namespace PacBio
