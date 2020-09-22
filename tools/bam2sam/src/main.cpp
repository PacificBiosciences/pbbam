// Author: Derek Barnett

#include <cstdlib>

#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include "Bam2SamSettings.h"
#include "Bam2SamWorkflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, PacBio::Bam2Sam::Settings::CreateCLI(),
                                   &PacBio::Bam2Sam::Workflow::Runner);
    } catch (const std::exception& e) {
        std::cerr << "bam2sam ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
