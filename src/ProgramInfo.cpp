#include "PbbamInternalConfig.h"

#include <pbbam/ProgramInfo.h>

#include <sstream>
#include <type_traits>

#include <pbbam/SamTagCodec.h>
#include <pbbam/StringUtilities.h>

namespace PacBio {
namespace BAM {
namespace {

const std::string ProgramInfoTokenID{"ID"};
const std::string ProgramInfoTokenCL{"CL"};
const std::string ProgramInfoTokenDS{"DS"};
const std::string ProgramInfoTokenPN{"PN"};
const std::string ProgramInfoTokenPP{"PP"};
const std::string ProgramInfoTokenVN{"VN"};

}  // namespace

ProgramInfo::ProgramInfo(std::string id) : id_{std::move(id)} {}

std::string ProgramInfo::CommandLine() const { return commandLine_; }

ProgramInfo& ProgramInfo::CommandLine(std::string cmd)
{
    commandLine_ = std::move(cmd);
    return *this;
}

std::map<std::string, std::string> ProgramInfo::CustomTags() const { return custom_; }

ProgramInfo& ProgramInfo::CustomTags(std::map<std::string, std::string> custom)
{
    custom_ = std::move(custom);
    return *this;
}

std::string ProgramInfo::Description() const { return description_; }

ProgramInfo& ProgramInfo::Description(std::string description)
{
    description_ = std::move(description);
    return *this;
}

ProgramInfo ProgramInfo::FromSam(const std::string& sam)
{
    // pop off '@PG\t', then split rest of line into tokens
    const auto tokens = Split(sam.substr(4), '\t');
    if (tokens.empty()) {
        return {};
    }

    ProgramInfo prog;
    std::map<std::string, std::string> custom;

    // iterate over tokens
    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set program contents
        if (tokenTag == ProgramInfoTokenID) {
            prog.Id(std::move(tokenValue));
        } else if (tokenTag == ProgramInfoTokenCL) {
            prog.CommandLine(std::move(tokenValue));
        } else if (tokenTag == ProgramInfoTokenDS) {
            prog.Description(std::move(tokenValue));
        } else if (tokenTag == ProgramInfoTokenPN) {
            prog.Name(std::move(tokenValue));
        } else if (tokenTag == ProgramInfoTokenPP) {
            prog.PreviousProgramId(std::move(tokenValue));
        } else if (tokenTag == ProgramInfoTokenVN) {
            prog.Version(std::move(tokenValue));

            // otherwise, "custom" tag
        } else {
            custom[tokenTag] = std::move(tokenValue);
        }
    }

    prog.CustomTags(custom);
    return prog;
}

std::string ProgramInfo::Id() const { return id_; }

ProgramInfo& ProgramInfo::Id(std::string id)
{
    id_ = std::move(id);
    return *this;
}

bool ProgramInfo::IsValid() const { return !id_.empty(); }

std::string ProgramInfo::Name() const { return name_; }

ProgramInfo& ProgramInfo::Name(std::string name)
{
    name_ = std::move(name);
    return *this;
}

std::string ProgramInfo::PreviousProgramId() const { return previousProgramId_; }

ProgramInfo& ProgramInfo::PreviousProgramId(std::string id)
{
    previousProgramId_ = std::move(id);
    return *this;
}

std::string ProgramInfo::ToSam(const ProgramInfo& prog) { return prog.ToSam(); }

std::string ProgramInfo::ToSam() const
{
    std::ostringstream out;
    out << "@PG" << MakeSamTag(ProgramInfoTokenID, id_);

    if (!name_.empty()) {
        out << MakeSamTag(ProgramInfoTokenPN, name_);
    }
    if (!version_.empty()) {
        out << MakeSamTag(ProgramInfoTokenVN, version_);
    }
    if (!description_.empty()) {
        out << MakeSamTag(ProgramInfoTokenDS, description_);
    }
    if (!previousProgramId_.empty()) {
        out << MakeSamTag(ProgramInfoTokenPP, previousProgramId_);
    }
    if (!commandLine_.empty()) {
        out << MakeSamTag(ProgramInfoTokenCL, commandLine_);
    }

    // append any custom tags
    for (const auto& attribute : custom_) {
        out << MakeSamTag(attribute.first, attribute.second);
    }
    return out.str();
}

std::string ProgramInfo::Version() const { return version_; }

ProgramInfo& ProgramInfo::Version(std::string version)
{
    version_ = std::move(version);
    return *this;
}

}  // namespace BAM
}  // namespace PacBio
