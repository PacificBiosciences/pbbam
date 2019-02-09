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

#include "pbbam/SamTagCodec.h"
#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {
namespace {

const std::string token_SN{"SN"};
const std::string token_LN{"LN"};
const std::string token_AS{"AS"};
const std::string token_M5{"M5"};
const std::string token_SP{"SP"};
const std::string token_UR{"UR"};

}  // anonymous

SequenceInfo::SequenceInfo(std::string name, std::string length)
    : name_(std::move(name)), length_(std::move(length))
{
}

bool SequenceInfo::operator==(const SequenceInfo& other) const
{
    return assemblyId_ == other.assemblyId_ && checksum_ == other.checksum_ &&
           length_ == other.length_ && name_ == other.name_ && species_ == other.species_ &&
           uri_ == other.uri_ && custom_ == other.custom_;
}

bool SequenceInfo::operator!=(const SequenceInfo& other) const { return !(*this == other); }

std::string SequenceInfo::AssemblyId() const { return assemblyId_; }

SequenceInfo& SequenceInfo::AssemblyId(std::string id)
{
    assemblyId_ = std::move(id);
    return *this;
}

std::string SequenceInfo::Checksum() const { return checksum_; }

SequenceInfo& SequenceInfo::Checksum(std::string checksum)
{
    checksum_ = std::move(checksum);
    return *this;
}

std::map<std::string, std::string> SequenceInfo::CustomTags() const { return custom_; }

SequenceInfo& SequenceInfo::CustomTags(std::map<std::string, std::string> custom)
{
    custom_ = std::move(custom);
    return *this;
}

SequenceInfo SequenceInfo::FromSam(const std::string& sam)
{
    // pop off '@SQ\t', then split rest of line into tokens
    const auto tokens = Split(sam.substr(4), '\t');
    if (tokens.empty()) return {};

    SequenceInfo seq;
    std::map<std::string, std::string> custom;

    // iterate over tokens
    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set sequence info
        // clang-format off
        if      (tokenTag == token_SN) seq.Name(std::move(tokenValue));
        else if (tokenTag == token_LN) seq.Length(std::move(tokenValue));
        else if (tokenTag == token_AS) seq.AssemblyId(std::move(tokenValue));
        else if (tokenTag == token_M5) seq.Checksum(std::move(tokenValue));
        else if (tokenTag == token_SP) seq.Species(std::move(tokenValue));
        else if (tokenTag == token_UR) seq.Uri(std::move(tokenValue));
        // clang-format on

        // otherwise, "custom" tag
        else
            custom[tokenTag] = std::move(tokenValue);
    }

    seq.CustomTags(std::move(custom));
    return seq;
}

bool SequenceInfo::IsValid() const
{
    if (name_.empty()) return false;

    // use long instead of int32_t, just to make sure we can catch overflow
    const long l = atol(length_.c_str());
    return l >= 0 && l <= std::numeric_limits<int32_t>::max();
}

std::string SequenceInfo::Length() const { return length_; }

SequenceInfo& SequenceInfo::Length(std::string length)
{
    length_ = std::move(length);
    return *this;
}

std::string SequenceInfo::Name() const { return name_; }

SequenceInfo& SequenceInfo::Name(std::string name)
{
    name_ = std::move(name);
    return *this;
}

std::string SequenceInfo::Species() const { return species_; }

SequenceInfo& SequenceInfo::Species(std::string species)
{
    species_ = std::move(species);
    return *this;
}

std::string SequenceInfo::ToSam(const SequenceInfo& seq) { return seq.ToSam(); }

std::string SequenceInfo::ToSam() const
{
    std::ostringstream out;
    out << "@SQ" << MakeSamTag(token_SN, name_);

    // clang-format off
    if (!length_.empty())     out << MakeSamTag(token_LN, length_);
    if (!assemblyId_.empty()) out << MakeSamTag(token_AS, assemblyId_);
    if (!checksum_.empty())   out << MakeSamTag(token_M5, checksum_);
    if (!species_.empty())    out << MakeSamTag(token_SP, species_);
    if (!uri_.empty())        out << MakeSamTag(token_UR, uri_);
    // clang-format on

    // append any custom tags
    for (auto&& attribute : custom_)
        out << MakeSamTag(std::move(attribute.first), std::move(attribute.second));

    return out.str();
}

std::string SequenceInfo::Uri() const { return uri_; }

SequenceInfo& SequenceInfo::Uri(std::string uri)
{
    uri_ = std::move(uri);
    return *this;
}

}  // namespace BAM
}  // namespace PacBio
