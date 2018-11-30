// File Description
/// \file ProgramInfo.cpp
/// \brief Implements the ProgramInfo class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ProgramInfo.h"

#include <sstream>

#include "pbbam/SamTagCodec.h"
#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {
namespace {

const std::string ProgramInfoTokenID{"ID"};
const std::string ProgramInfoTokenCL{"CL"};
const std::string ProgramInfoTokenDS{"DS"};
const std::string ProgramInfoTokenPN{"PN"};
const std::string ProgramInfoTokenPP{"PP"};
const std::string ProgramInfoTokenVN{"VN"};

}  // anonymous

ProgramInfo::ProgramInfo(std::string id) : id_{std::move(id)} {}

ProgramInfo ProgramInfo::FromSam(const std::string& sam)
{
    // pop off '@PG\t', then split rest of line into tokens
    const auto tokens = Split(sam.substr(4), '\t');
    if (tokens.empty()) return {};

    ProgramInfo prog;
    std::map<std::string, std::string> custom;

    // iterate over tokens
    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set program contents
        // clang-format off
        if      (tokenTag == ProgramInfoTokenID) prog.Id(std::move(tokenValue));
        else if (tokenTag == ProgramInfoTokenCL) prog.CommandLine(std::move(tokenValue));
        else if (tokenTag == ProgramInfoTokenDS) prog.Description(std::move(tokenValue));
        else if (tokenTag == ProgramInfoTokenPN) prog.Name(std::move(tokenValue));
        else if (tokenTag == ProgramInfoTokenPP) prog.PreviousProgramId(std::move(tokenValue));
        else if (tokenTag == ProgramInfoTokenVN) prog.Version(std::move(tokenValue));
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
    out << "@PG" << MakeSamTag(ProgramInfoTokenID, id_);

    // clang-format off
    if (!name_.empty())              out << MakeSamTag(ProgramInfoTokenPN, name_);
    if (!version_.empty())           out << MakeSamTag(ProgramInfoTokenVN, version_);
    if (!description_.empty())       out << MakeSamTag(ProgramInfoTokenDS, description_);
    if (!previousProgramId_.empty()) out << MakeSamTag(ProgramInfoTokenPP, previousProgramId_);
    if (!commandLine_.empty())       out << MakeSamTag(ProgramInfoTokenCL, commandLine_);
    // clang-format on

    // append any custom tags
    for (const auto& attribute : custom_)
        out << MakeSamTag(attribute.first, attribute.second);
    return out.str();
}

}  // namespace BAM
}  // namespace PacBio
