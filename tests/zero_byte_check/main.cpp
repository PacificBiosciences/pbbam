#include <iostream>
#include <stdexcept>
#include <tuple>

#include <pbbam/BamReader.h>
#include <pbbam/SamReader.h>

int main(int /* argc */, char* argv[])
{
    try {
        if (std::string{argv[1]} == "bam") {
            PacBio::BAM::BamReader reader{"-"};
            std::ignore = reader;
        } else {
            PacBio::BAM::SamReader reader{"-"};
            std::ignore = reader;
        }
    } catch (const std::exception& e) {
        std::cout << e.what();
    }
    return 0;
}
