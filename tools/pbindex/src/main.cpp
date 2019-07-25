// Author: Derek Barnett

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include "PbIndexSettings.h"
#include "PbIndexWorkflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, PacBio::PbIndex::Settings::CreateCLI(),
                                   &PacBio::PbIndex::Workflow::Runner);
    } catch (const std::exception& e) {
        std::cerr << "pbindex ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
