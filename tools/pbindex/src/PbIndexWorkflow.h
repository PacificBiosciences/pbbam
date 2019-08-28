// Author: Derek Barnett

#ifndef PBINDEX_WORKFLOW_H
#define PBINDEX_WORKFLOW_H

#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace PbIndex {

struct Workflow
{
    static int Runner(const CLI_v2::Results& args);
};

}  // namespace PbIndex
}  // namespace PacBio

#endif  // PBINDEX_WORKFLOW_H
