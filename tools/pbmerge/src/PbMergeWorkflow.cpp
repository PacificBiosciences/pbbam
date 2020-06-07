// Author: Derek Barnett

#include "PbMergeWorkflow.h"

#include <pbbam/BamFileMerger.h>
#include <pbbam/DataSet.h>
#include <pbbam/ProgramInfo.h>

#include "PbMergeSettings.h"
#include "PbMergeVersion.h"

namespace PacBio {
namespace PbMerge {

int Workflow::Runner(const CLI_v2::Results& args)
{
    const Settings settings{args};

    // setup our @PG entry to add to header
    BAM::ProgramInfo mergeProgram;
    mergeProgram.Id(std::string("pbmerge-") + PbMerge::Version)
        .Name("pbmerge")
        .Version(PbMerge::Version);

    BAM::DataSet dataset;
    if (settings.InputFiles.size() == 1)
        dataset = BAM::DataSet(settings.InputFiles.front());
    else
        dataset = BAM::DataSet(settings.InputFiles);

    BAM::BamFileMerger::Merge(dataset, settings.OutputFile, settings.CreatePbi, mergeProgram);

    return EXIT_SUCCESS;
}

}  // namespace PbMerge
}  // namespace PacBio
