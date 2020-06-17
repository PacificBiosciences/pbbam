// File Description
/// \file BamRecordTags.h
/// \brief Defines the BamRecordTags utility class.
//
// Author: Derek Barnett

#ifndef BAMRECORDTAGS_H
#define BAMRECORDTAGS_H

#include "pbbam/Config.h"

#include <cassert>

#include <string>
#include <unordered_map>

#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordImpl.h>
#include <pbbam/BamRecordTag.h>

namespace PacBio {
namespace BAM {

class BamRecordTags
{
public:
    // tag info
    static bool IsPulse(const BamRecordTag tag);
    static std::string LabelFor(const BamRecordTag tag);

private:
    struct BamRecordTagData
    {
        const std::string label_;  //[3]; // 2-char tag plus NULL
        const bool isPulse_;
    };

    using TagLookupType = std::unordered_map<BamRecordTag, BamRecordTagData>;
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

}  // namespace BAM
}  // namespace PacBio

#endif  // BAMRECORDTAGS_H
