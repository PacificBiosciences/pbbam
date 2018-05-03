// File Description
/// \file ProgramInfo.cpp
/// \brief Implements the ProgramInfo class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ProgramInfo.h"

#include <sstream>

#include "SequenceUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

static std::string ProgramInfoTokenID{"ID"};
static std::string ProgramInfoTokenCL{"CL"};
static std::string ProgramInfoTokenDS{"DS"};
static std::string ProgramInfoTokenPN{"PN"};
static std::string ProgramInfoTokenPP{"PP"};
static std::string ProgramInfoTokenVN{"VN"};

}  // namespace internal

ProgramInfo::ProgramInfo(std::string id) : id_{std::move(id)} {}

ProgramInfo ProgramInfo::FromSam(const std::string& sam)
{
    // pop off '@PG\t', then split rest of line into tokens
    const auto tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty()) return {};

    ProgramInfo prog;
    std::map<std::string, std::string> custom;

    // iterate over tokens
    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set program contents
        // clang-format off
        if      (tokenTag == internal::ProgramInfoTokenID) prog.Id(std::move(tokenValue));
        else if (tokenTag == internal::ProgramInfoTokenCL) prog.CommandLine(std::move(tokenValue));
        else if (tokenTag == internal::ProgramInfoTokenDS) prog.Description(std::move(tokenValue));
        else if (tokenTag == internal::ProgramInfoTokenPN) prog.Name(std::move(tokenValue));
        else if (tokenTag == internal::ProgramInfoTokenPP) prog.PreviousProgramId(std::move(tokenValue));
        else if (tokenTag == internal::ProgramInfoTokenVN) prog.Version(std::move(tokenValue));
        // clang-format on

        // otherwise, "custom" tag
        else
            custom[tokenTag] = std::move(tokenValue);
    }

    prog.CustomTags(custom);
    return prog;
}

std::string ProgramInfo::ToSam() const
{
    std::ostringstream out;
    out << "@PG" << internal::MakeSamTag(internal::ProgramInfoTokenID, id_);

    // clang-format off
    if (!name_.empty())              out << internal::MakeSamTag(internal::ProgramInfoTokenPN, name_);
    if (!version_.empty())           out << internal::MakeSamTag(internal::ProgramInfoTokenVN, version_);
    if (!description_.empty())       out << internal::MakeSamTag(internal::ProgramInfoTokenDS, description_);
    if (!previousProgramId_.empty()) out << internal::MakeSamTag(internal::ProgramInfoTokenPP, previousProgramId_);
    if (!commandLine_.empty())       out << internal::MakeSamTag(internal::ProgramInfoTokenCL, commandLine_);
    // clang-format on

    // append any custom tags
    for (const auto& attribute : custom_)
        out << internal::MakeSamTag(attribute.first, attribute.second);
    return out.str();
}

}  // namespace BAM
}  // namespace PacBio
