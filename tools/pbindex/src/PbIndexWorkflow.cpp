// Author: Derek Barnett

#include "PbIndexWorkflow.h"

#include <pbbam/BamFile.h>
#include <pbbam/PbiFile.h>

#include "PbIndexSettings.h"

namespace PacBio {
namespace PbIndex {

int Workflow::Runner(const CLI_v2::Results& args)
{
    const Settings settings{args};
    BAM::PbiFile::CreateFrom(settings.InputFile);
    return EXIT_SUCCESS;
}

}  // namespace PbIndex
}  // namespace PacBio