#ifndef CCSKINETICSBYSTRANDIFY_WORKFLOW_H
#define CCSKINETICSBYSTRANDIFY_WORKFLOW_H

#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace CcsKineticsBystrandify {

struct Workflow
{
    ///
    /// \brief Takes a PacBio CCS BAM and converts it into a pseudo-subreads
    ///        bystrand CCS look-a-like BAM.
    ///
    /// \throws std::runtime_error
    ///
    static int Runner(const CLI_v2::Results& args);
};

}  // namespace CcsKineticsBystrandify
}  // namespace PacBio

#endif  // CCSKINETICSBYSTRANDIFY_WORKFLOW_H
