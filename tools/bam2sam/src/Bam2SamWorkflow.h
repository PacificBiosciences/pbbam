// Author: Derek Barnett

#ifndef BAM2SAM_WORKFLOW_H
#define BAM2SAM_WORKFLOW_H

#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Bam2Sam {

struct Workflow
{
    static int Runner(const CLI_v2::Results& args);
};

}  // namespace Bam2Sam
}  // namespace PacBio

#endif  // BAM2SAM_WORKFLOW_H
