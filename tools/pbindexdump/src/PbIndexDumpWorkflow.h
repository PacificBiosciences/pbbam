// Author: Derek Barnett

#ifndef PBINDEXDUMP_WORKFLOW_H
#define PBINDEXDUMP_WORKFLOW_H

#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace PbIndexDump {

struct Workflow
{
    static int Runner(const CLI_v2::Results& args);
};

}  // namespace PbIndexDump
}  // namespace PacBio

#endif  // PBINDEXDUMP_WORKFLOW_H
