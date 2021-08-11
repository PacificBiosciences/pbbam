#include <pbbam/../../src/TimeUtils.h>

#include <ctime>

#include <chrono>
#include <string>

#include <gtest/gtest.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_TimeUtils, can_convert_time_point_to_iso8601)
{
    const time_t rawTime = 436428750L;
    const auto timestamp = std::chrono::system_clock::from_time_t(rawTime);

    const std::string expected = "1983-10-31T06:12:30Z";  // no ms in test case
    const auto actual = TimeUtils::ToIso8601(timestamp);
    EXPECT_EQ(expected, actual);
}

TEST(BAM_TimeUtils, can_convert_time_point_to_pacbio_dataset_format)
{
    const time_t rawTime = 436428750L;
    const auto timestamp = std::chrono::system_clock::from_time_t(rawTime);

    const std::string expected = std::string{"831031_061230"};  // no ms in test case
    const std::string actual = TimeUtils::ToDataSetFormat(timestamp);
    EXPECT_EQ(expected, actual);
}
