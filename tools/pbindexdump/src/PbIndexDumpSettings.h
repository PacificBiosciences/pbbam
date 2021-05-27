#ifndef PBINDEXDUMP_SETTINGS_H
#define PBINDEXDUMP_SETTINGS_H

#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace PbIndexDump {

struct Settings
{
    static CLI_v2::Interface CreateCLI();

    explicit Settings(const CLI_v2::Results& args);

    std::string InputFile;
    std::string Format;
    int JsonIndentLevel = 4;
    bool JsonRaw = false;
};

}  // namespace PbIndexDump
}  // namespace PacBio

#endif  // PBINDEXDUMP_SETTINGS_H
