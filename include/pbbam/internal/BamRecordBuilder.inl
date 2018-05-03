// File Description
/// \file BamRecordBuilder.inl
/// \brief Inline implementations for the BamRecordBuilder class.
//
// Author: Derek Barnett

#include "pbbam/BamRecordBuilder.h"

namespace PacBio {
namespace BAM {

inline BamRecordBuilder& BamRecordBuilder::Bin(const uint32_t bin)
{ core_.bin = bin; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Flag(const uint32_t flag)
{ core_.flag = flag; return *this; }

inline BamRecordBuilder& BamRecordBuilder::InsertSize(const int32_t iSize)
{ core_.isize = iSize; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MapQuality(const uint8_t mapQual)
{ core_.qual = mapQual; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MatePosition(const int32_t pos)
{ core_.mpos = pos; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MateReferenceId(const int32_t id)
{ core_.mtid = id; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Position(const int32_t pos)
{ core_.pos = pos; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Qualities(std::string qualities)
{ qualities_ = std::move(qualities); return *this; }

inline BamRecordBuilder& BamRecordBuilder::ReferenceId(const int32_t id)
{ core_.tid = id; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Tags(TagCollection tags)
{ tags_ = std::move(tags); return *this; }

} // namespace BAM
} // namespace PacBio
