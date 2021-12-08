#ifndef PBBAM_BAMRECORDTAGS_H
#define PBBAM_BAMRECORDTAGS_H

#include <pbbam/Config.h>

#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordImpl.h>
#include <pbbam/BamRecordTag.h>

#include <string>
#include <unordered_map>

#include <cassert>

namespace PacBio {
namespace BAM {

class BamRecordTags
{
public:
    // tag info
    static bool IsIPD(const BamRecordTag tag);
    static bool IsPulse(const BamRecordTag tag);
    static bool IsPW(const BamRecordTag tag);
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

inline bool BamRecordTags::IsIPD(const BamRecordTag tag)
{
    return tag == BamRecordTag::IPD || tag == BamRecordTag::FORWARD_IPD ||
           tag == BamRecordTag::REVERSE_IPD;
}

inline bool BamRecordTags::IsPulse(const BamRecordTag tag)
{
    assert(tagLookup.find(tag) != tagLookup.cend());
    return tagLookup.at(tag).isPulse_;
}

inline bool BamRecordTags::IsPW(const BamRecordTag tag)
{
    return tag == BamRecordTag::PULSE_WIDTH || tag == BamRecordTag::FORWARD_PW ||
           tag == BamRecordTag::REVERSE_PW;
}

inline std::string BamRecordTags::LabelFor(const BamRecordTag tag)
{
    assert(tagLookup.find(tag) != tagLookup.cend());
    return tagLookup.at(tag).label_;
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BAMRECORDTAGS_H
