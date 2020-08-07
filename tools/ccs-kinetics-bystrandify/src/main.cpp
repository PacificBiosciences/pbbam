#include <cstdlib>

#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include "CcsKineticsBystrandifySettings.h"
#include "CcsKineticsBystrandifyWorkflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv,
                                   PacBio::CcsKineticsBystrandify::Settings::CreateCLI(),
                                   &PacBio::CcsKineticsBystrandify::Workflow::Runner);
    } catch (const std::exception& e) {
        std::cerr << "[ccs-kinetics-bystrandify] ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
