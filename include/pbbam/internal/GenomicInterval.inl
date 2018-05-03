// File Description
/// \file GenomicInterval.inl
/// \brief Inline implementations for the GenomicInterval class.
//
// Author: Derek Barnett

#include "pbbam/GenomicInterval.h"

namespace PacBio {
namespace BAM {

inline std::string GenomicInterval::Name() const
{ return name_; }

inline GenomicInterval& GenomicInterval::Name(std::string name)
{ name_ = std::move(name); return *this; }

inline PacBio::BAM::Interval<Position> GenomicInterval::Interval() const
{ return interval_; }

inline GenomicInterval& GenomicInterval::Interval(PacBio::BAM::Interval<Position> interval)
{ interval_ = std::move(interval); return *this; }

inline bool GenomicInterval::IsValid() const
{
    return !name_.empty() &&
           interval_.Start() >= 0 &&
           interval_.Stop()  >= 0 &&
           interval_.IsValid();
}

inline size_t GenomicInterval::Length() const
{ return interval_.Length(); }

inline Position GenomicInterval::Start() const
{ return interval_.Start(); }

inline GenomicInterval& GenomicInterval::Start(const Position start)
{ interval_.Start(start); return *this; }

inline Position GenomicInterval::Stop() const
{ return interval_.Stop(); }

inline GenomicInterval& GenomicInterval::Stop(const Position stop)
{ interval_.Stop(stop); return *this; }

inline bool GenomicInterval::operator==(const GenomicInterval& other) const
{ return name_ == other.name_ && interval_ == other.interval_; }

inline bool GenomicInterval::operator!=(const GenomicInterval& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio
