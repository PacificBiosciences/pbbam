#ifndef CCSKINETICSBYSTRANDIFY_SETTINGS_H
#define CCSKINETICSBYSTRANDIFY_SETTINGS_H

#include <cstdint>

#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace CcsKineticsBystrandify {

struct Settings
{
    static CLI_v2::Interface CreateCLI();

    explicit Settings(const CLI_v2::Results& args);

    std::string InputFilename;
    std::string OutputFilename;

    int32_t MinCoverage;
};

}  // namespace CcsKineticsBystrandify
}  // namespace PacBio

#endif  // CCSKINETICSBYSTRANDIFY_SETTINGS_H
