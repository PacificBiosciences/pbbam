// Author: Derek Barnett

#include "PbIndex.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbbam/BamFile.h>
#include <pbbam/PbiRawData.h>

using namespace pbindex;

Settings::Settings() : printPbiContents_(false) {}

int PbIndex::Create(const Settings& settings)
{
    try {
        PacBio::BAM::BamFile bamFile(settings.inputBamFilename_);
        bamFile.CreatePacBioIndex();
        return EXIT_SUCCESS;
    } catch (std::runtime_error& e) {
        std::cerr << "pbindex ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}

//int PbIndex::Print(const Settings& settings)
//{

//}

int PbIndex::Run(const Settings& settings)
{
    //    if (settings.printPbiContents_)
    //        return Print(settings);
    //    else
    return Create(settings);
}
