#include "PbbamInternalConfig.h"

#include <pbbam/RecordType.h>

#include <map>
#include <stdexcept>

namespace PacBio {
namespace BAM {

bool IsCcsOrTranscript(const RecordType type)
{
    return (type == RecordType::CCS) || (type == RecordType::TRANSCRIPT);
}

RecordType RecordTypeFromString(const std::string& type)
{
    // clang-format off
    static const auto lookup = std::map<std::string, RecordType>
    {
        { "ZMW",        RecordType::ZMW },
        { "HQREGION",   RecordType::HQREGION },
        { "SUBREAD",    RecordType::SUBREAD },
        { "CCS",        RecordType::CCS },
        { "SCRAP",      RecordType::SCRAP },
        { "TRANSCRIPT", RecordType::TRANSCRIPT },
        { "UNKNOWN",    RecordType::UNKNOWN },
        { "SEGMENT",    RecordType::SEGMENT },
    };
    // clang-format on

    // match exising BamRecord behavior: fallback to UNKNOWN
    const auto found = lookup.find(type);
    if (found != lookup.cend()) {
        return found->second;
    }
    return RecordType::UNKNOWN;
}

std::string ToString(const RecordType type)
{
    // clang-format off
    static const auto lookup = std::map<RecordType, std::string>
    {
        { RecordType::ZMW,        "ZMW" },
        { RecordType::HQREGION,   "HQREGION" },
        { RecordType::SUBREAD,    "SUBREAD" },
        { RecordType::CCS,        "CCS" },
        { RecordType::SCRAP,      "SCRAP" },
        { RecordType::TRANSCRIPT, "TRANSCRIPT" },
        { RecordType::UNKNOWN,    "UNKNOWN" },
        { RecordType::SEGMENT,    "SEGMENT" },
    };
    // clang-format on

    try {
        return lookup.at(type);
    } catch (std::exception&) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: unknown record type"};
    }
}

}  // namespace BAM
}  // namespace PacBio
