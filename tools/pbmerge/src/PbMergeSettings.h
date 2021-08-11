#ifndef PBMERGE_SETTINGS_H
#define PBMERGE_SETTINGS_H

#include <string>
#include <vector>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace PbMerge {

struct Settings
{
    static CLI_v2::Interface CreateCLI();

    explicit Settings(const CLI_v2::Results& args);

    std::vector<std::string> InputFiles;
    std::string OutputFile;
    bool CreatePbi;
    std::vector<std::string> errors_;
};

}  // namespace PbMerge
}  // namespace PacBio

#endif  // PBMERGE_SETTINGS_H
