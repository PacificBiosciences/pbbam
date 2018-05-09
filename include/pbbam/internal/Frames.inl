// File Description
/// \file Frames.inl
/// \brief Inline implementations for the Frames class.
//
// Author: Derek Barnett

#include "pbbam/Frames.h"

namespace PacBio {
namespace BAM {

inline const std::vector<uint16_t>& Frames::Data() const
{ return data_; }

inline std::vector<uint16_t>& Frames::DataRaw()
{ return data_; }

inline std::vector<uint8_t> Frames::Encode() const
{ return Frames::Encode(data_); }

inline Frames& Frames::Data(std::vector<uint16_t> frames)
{ data_ = std::move(frames); return *this; }

inline std::vector<uint16_t>::const_iterator Frames::begin() const
{ return data_.begin(); }

inline std::vector<uint16_t>::iterator Frames::begin()
{ return data_.begin(); }

inline std::vector<uint16_t>::const_iterator Frames::cbegin() const
{ return data_.cbegin(); }

inline std::vector<uint16_t>::const_iterator Frames::cend() const
{ return data_.cend(); }

inline std::vector<uint16_t>::const_iterator Frames::end() const
{ return data_.end(); }

inline std::vector<uint16_t>::iterator Frames::end()
{ return data_.end(); }

inline size_t Frames::size() const
{ return data_.size(); }

inline bool Frames::empty() const
{ return data_.empty(); }

inline bool Frames::operator==(const Frames& other) const
{ return data_ == other.data_; }

inline bool Frames::operator!=(const Frames& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio
