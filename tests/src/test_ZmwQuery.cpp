// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include <pbbam/ZmwGroupQuery.h>
#include <pbbam/ZmwQuery.h>

#include "PbbamTestData.h"

namespace ZmwQueryTests {

const std::string input =
    PacBio::BAM::PbbamTestsConfig::Data_Dir + "/chunking/chunking.subreadset.xml";

}  // namespace ZmwQueryTests

TEST(ZmwQueryTest, whitelist_query_returns_nothing_from_empty_whitelist)
{
    const std::vector<int32_t> whitelist;

    size_t count = 0;
    PacBio::BAM::ZmwQuery query{whitelist, ZmwQueryTests::input};
    for (const auto& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(0, count);
}

TEST(ZmwQueryTest, whitelist_query_returns_only_requested_zmws)
{
    const std::vector<int32_t> whitelist{
        1411,   // 12
        54636,  // 26
        109697  // 10
    };

    size_t count = 0;
    PacBio::BAM::ZmwQuery query{whitelist, ZmwQueryTests::input};
    for (const auto& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(48, count);
}

TEST(ZmwGroupQueryTest, whitelist_query_returns_nothing_from_empty_whitelist)
{
    const std::vector<int32_t> whitelist;

    size_t zmwCount = 0;
    size_t recordCount = 0;
    PacBio::BAM::ZmwGroupQuery query{whitelist, ZmwQueryTests::input};
    for (const auto& zmw : query) {
        ++zmwCount;
        for (const auto& record : zmw) {
            (void)record;
            ++recordCount;
        }
    }
    EXPECT_EQ(0, zmwCount);
    EXPECT_EQ(0, recordCount);
}

TEST(ZmwGroupQueryTest, whitelist_query_returns_only_requested_zmws)
{
    const std::vector<int32_t> whitelist{
        1411,   // 12
        54636,  // 26
        109697  // 10
    };

    size_t zmwCount = 0;
    size_t recordCount = 0;
    PacBio::BAM::ZmwGroupQuery query{whitelist, ZmwQueryTests::input};
    for (const auto& zmw : query) {
        ++zmwCount;
        for (const auto& record : zmw) {
            (void)record;
            ++recordCount;
        }
    }
    EXPECT_EQ(3, zmwCount);
    EXPECT_EQ(48, recordCount);
}

TEST(ZmwGroupQueryTest, round_robin_query_can_return_records_applying_dataset_filter)
{
    size_t zmwCount = 0;
    size_t recordCount = 0;
    PacBio::BAM::ZmwGroupQuery query{ZmwQueryTests::input, PacBio::BAM::DataSetFilterMode::APPLY};
    for (const auto& zmw : query) {
        ++zmwCount;
        for (const auto& record : zmw) {
            (void)record;
            ++recordCount;
        }
    }

    // zmw < 1816
    EXPECT_EQ(15, zmwCount);
    EXPECT_EQ(150, recordCount);
}

TEST(ZmwGroupQueryTest, round_robin_query_can_return_records_ignoring_dataset_filter)
{
    size_t zmwCount = 0;
    size_t recordCount = 0;
    PacBio::BAM::ZmwGroupQuery query{ZmwQueryTests::input, PacBio::BAM::DataSetFilterMode::IGNORE};
    for (const auto& zmw : query) {
        ++zmwCount;
        for (const auto& record : zmw) {
            (void)record;
            ++recordCount;
        }
    }
    EXPECT_EQ(90, zmwCount);       // 30 + 30 + 30
    EXPECT_EQ(1220, recordCount);  // 432 + 409 + 379
}
