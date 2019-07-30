// Author: Derek Barnett

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include "PbMergeSettings.h"
#include "PbMergeWorkflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, PacBio::PbMerge::Settings::CreateCLI(),
                                   &PacBio::PbMerge::Workflow::Runner);
    } catch (const std::exception& e) {
        std::cerr << "pbmerge ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
