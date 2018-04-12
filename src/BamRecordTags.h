// File Description
/// \file BamRecordTags.h
/// \brief Defines the BamRecordTags utility class.
//
// Author: Derek Barnett

#ifndef BAMRECORDTAGS_H
#define BAMRECORDTAGS_H

#include <cassert>
#include <string>
#include <unordered_map>

#include "EnumClassHash.h"
#include "pbbam/BamRecord.h"
#include "pbbam/BamRecordImpl.h"
#include "pbbam/BamRecordTag.h"

namespace PacBio {
namespace BAM {
namespace internal {

class BamRecordTags
{
public:
    // tag info
    static inline bool IsPulse(const BamRecordTag tag);
    static inline std::string LabelFor(const BamRecordTag tag);

private:
    struct BamRecordTagData
    {
        const std::string label_;  //[3]; // 2-char tag plus NULL
        const bool isPulse_;
    };
    typedef std::unordered_map<BamRecordTag, BamRecordTagData, EnumClassHash> TagLookupType;

    static const TagLookupType tagLookup;
};

inline bool BamRecordTags::IsPulse(const BamRecordTag tag)
{
    assert(tagLookup.find(tag) != tagLookup.cend());
    return tagLookup.at(tag).isPulse_;
}

inline std::string BamRecordTags::LabelFor(const BamRecordTag tag)
{
    assert(tagLookup.find(tag) != tagLookup.cend());
    return tagLookup.at(tag).label_;
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // BAMRECORDTAGS_H
