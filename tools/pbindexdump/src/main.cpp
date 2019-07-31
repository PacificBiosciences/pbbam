// Author: Derek Barnett

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include "PbIndexDumpSettings.h"
#include "PbIndexDumpWorkflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, PacBio::PbIndexDump::Settings::CreateCLI(),
                                   &PacBio::PbIndexDump::Workflow::Runner);
    } catch (const std::exception& e) {
        std::cerr << "pbindexdump ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
