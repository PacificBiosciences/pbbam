// Author: Derek Barnett

#ifndef PBMERGE_WORKFLOW_H
#define PBMERGE_WORKFLOW_H

#include <string>
#include <vector>

#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace PbMerge {

struct Workflow
{
    static int Runner(const CLI_v2::Results& args);
};

}  // namespace PbMerge
}  // namespace PacBio

#endif  // PBMERGE_WORKFLOW_H
