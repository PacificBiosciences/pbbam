// Author: Ivan Sovic

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include "PbBamifySettings.h"
#include "PbBamifyWorkflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, PacBio::PbBamify::Settings::CreateCLI(),
                                   &PacBio::PbBamify::Workflow::Runner);
    } catch (const std::exception& e) {
        std::cerr << "pbbamify ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
