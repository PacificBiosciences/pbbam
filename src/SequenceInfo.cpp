// File Description
/// \file SequenceInfo.cpp
/// \brief Implements the SequenceInfo class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SequenceInfo.h"

#include <cstdint>
#include <limits>
#include <sstream>

#include "SequenceUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

static const std::string token_SN{"SN"};
static const std::string token_LN{"LN"};
static const std::string token_AS{"AS"};
static const std::string token_M5 = std::string("M5");
static const std::string token_SP = std::string("SP");
static const std::string token_UR = std::string("UR");

}  // namespace internal

SequenceInfo::SequenceInfo(std::string name, std::string length)
    : name_(std::move(name)), length_(std::move(length))
{
}

SequenceInfo SequenceInfo::FromSam(const std::string& sam)
{
    // pop off '@SQ\t', then split rest of line into tokens
    const auto tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty()) return SequenceInfo();

    SequenceInfo seq;
    std::map<std::string, std::string> custom;

    // iterate over tokens
    for (const std::string& token : tokens) {
        const std::string& tokenTag = token.substr(0, 2);
        const std::string& tokenValue = token.substr(3);

        // set sequence info
        // clang-format off
        if      (tokenTag == internal::token_SN) seq.Name(tokenValue);
        else if (tokenTag == internal::token_LN) seq.Length(tokenValue);
        else if (tokenTag == internal::token_AS) seq.AssemblyId(tokenValue);
        else if (tokenTag == internal::token_M5) seq.Checksum(tokenValue);
        else if (tokenTag == internal::token_SP) seq.Species(tokenValue);
        else if (tokenTag == internal::token_UR) seq.Uri(tokenValue);
        // clang-format on

        // otherwise, "custom" tag
        else
            custom[tokenTag] = tokenValue;
    }

    seq.CustomTags(custom);
    return seq;
}

bool SequenceInfo::IsValid() const
{
    if (name_.empty()) return false;

    // use long instead of int32_t, just to make sure we can catch overflow
    const long l = atol(length_.c_str());
    return l >= 0 && l <= std::numeric_limits<int32_t>::max();
}

std::string SequenceInfo::ToSam() const
{
    std::stringstream out;
    out << "@SQ" << internal::MakeSamTag(internal::token_SN, name_);

    // clang-format off
    if (!length_.empty())     out << internal::MakeSamTag(internal::token_LN, length_);
    if (!assemblyId_.empty()) out << internal::MakeSamTag(internal::token_AS, assemblyId_);
    if (!checksum_.empty())   out << internal::MakeSamTag(internal::token_M5, checksum_);
    if (!species_.empty())    out << internal::MakeSamTag(internal::token_SP, species_);
    if (!uri_.empty())        out << internal::MakeSamTag(internal::token_UR, uri_);
    // clang-format on

    // append any custom tags
    for (const auto& attribute : custom_)
        out << internal::MakeSamTag(attribute.first, attribute.second);

    return out.str();
}

}  // namespace BAM
}  // namespace PacBio
