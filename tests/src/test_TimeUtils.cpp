// Author: Derek Barnett

#include <chrono>
#include <ctime>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/../../src/TimeUtils.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;

TEST(TimeUtilsTest, ToIso8601)
{
    const time_t rawTime = 436428750L;
    const auto timestamp = std::chrono::system_clock::from_time_t(rawTime);

    const auto expected = std::string{"1983-10-31T06:12:30Z"};  // no ms in test case
    const auto actual = internal::ToIso8601(timestamp);
    EXPECT_EQ(expected, actual);
}

TEST(TimeUtilsTest, ToDataSetFormat)
{
    const time_t rawTime = 436428750L;
    const auto timestamp = std::chrono::system_clock::from_time_t(rawTime);

    const auto expected = std::string{"831031_061230"};  // no ms in test case
    const std::string actual = internal::ToDataSetFormat(timestamp);
    EXPECT_EQ(expected, actual);
}
