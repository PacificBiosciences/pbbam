#ifndef BAM2SAM_SETTINGS_H
#define BAM2SAM_SETTINGS_H

#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Bam2Sam {

struct Settings
{
    static CLI_v2::Interface CreateCLI();

    explicit Settings(const CLI_v2::Results& args);

    std::string InputFilename;
    bool NoHeader = false;
    bool HeaderOnly = false;
};

}  // namespace Bam2Sam
}  // namespace PacBio

#endif  // BAM2SAM_SETTINGS_H
