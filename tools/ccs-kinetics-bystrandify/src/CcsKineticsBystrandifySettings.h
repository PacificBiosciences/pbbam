#ifndef CCSKINETICSBYSTRANDIFY_SETTINGS_H
#define CCSKINETICSBYSTRANDIFY_SETTINGS_H

#include <cstdint>

#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace CcsKineticsBystrandify {

struct Settings
{
    struct Defaults
    {
        static constexpr int32_t MinCoverage = 1;
    };

    std::string CLI;
    std::string InputFilename;
    std::string OutputFilename;

    int32_t MinCoverage = Defaults::MinCoverage;

    static CLI_v2::Interface CreateCLI();
    explicit Settings(const CLI_v2::Results& args);
};

}  // namespace CcsKineticsBystrandify
}  // namespace PacBio

#endif  // CCSKINETICSBYSTRANDIFY_SETTINGS_H
