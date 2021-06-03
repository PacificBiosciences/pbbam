#ifndef PBBAMIFY_SETTINGS_H
#define PBBAMIFY_SETTINGS_H

#include <cstdint>
#include <string>
#include <vector>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace PbBamify {

struct Settings
{
    static CLI_v2::Interface CreateCLI();

    explicit Settings(const CLI_v2::Results& args);

    std::string InputFilename;
    std::string OutputFilename;
    std::string ReferenceFilename;
    std::string PbbamFilename;
    std::vector<std::string> Errors;
    int32_t VerboseLevel;
};

}  // namespace PbBamify
}  // namespace PacBio

#endif  // PBBAMIFY_SETTINGS_H
