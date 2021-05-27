#include <pbbam/BamReader.h>

#include <sstream>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;

TEST(BAM_BamReader, handles_zero_byte_file)
{
    try {
        BAM::BamReader reader{BAM::PbbamTestsConfig::Data_Dir + "/zero_bytes.bam"};
        ASSERT_FALSE("should not get here");
    } catch (const std::exception& e) {
        const std::string msg{e.what()};
        EXPECT_TRUE(msg.find("[pbbam] BAM reader ERROR: could not read from empty input:") !=
                    std::string::npos);
    }
}
