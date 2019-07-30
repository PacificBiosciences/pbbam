// Author: Derek Barnett

#ifndef PBINDEX_SETTINGS_H
#define PBINDEX_SETTINGS_H

#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace PbIndex {

struct Settings
{
    static CLI_v2::Interface CreateCLI();

    explicit Settings(const CLI_v2::Results& args);

    std::string InputFile;
};

}  // namespace PbIndex
}  // namespace PacBio

#endif  // PBINDEX_SETTINGS_H
