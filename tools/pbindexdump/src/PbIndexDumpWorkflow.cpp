#include "PbIndexDumpWorkflow.h"

#include <cassert>
#include <cstdlib>

#include "CppFormatter.h"
#include "JsonFormatter.h"
#include "PbIndexDumpSettings.h"

namespace PacBio {
namespace PbIndexDump {

int Workflow::Runner(const CLI_v2::Results& args)
{
    const Settings settings{args};
    if (settings.Format == "json") {
        JsonFormatter::Run(settings);
    } else {
        assert(settings.Format == "cpp");
        CppFormatter::Run(settings);
    }
    return EXIT_SUCCESS;
}

}  // namespace PbIndexDump
}  // namespace PacBio
