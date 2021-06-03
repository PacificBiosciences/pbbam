#include <pbbam/ZmwQuery.h>

#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/ZmwGroupQuery.h>

#include "PbbamTestData.h"

namespace ZmwQueryTests {

const std::string input =
    PacBio::BAM::PbbamTestsConfig::Data_Dir + "/chunking/chunking.subreadset.xml";

}  // namespace ZmwQueryTests

TEST(BAM_ZmwQuery, whitelist_query_returns_nothing_from_empty_whitelist)
{
    const std::vector<int32_t> whitelist;

    PacBio::BAM::ZmwQuery query{whitelist, ZmwQueryTests::input};
    EXPECT_EQ(0, std::distance(query.begin(), query.end()));
}

TEST(BAM_ZmwQuery, whitelist_query_returns_only_requested_zmws)
{
    const std::vector<int32_t> whitelist{
        1411,   // 12
        54636,  // 26
        109697  // 10
    };

    PacBio::BAM::ZmwQuery query{whitelist, ZmwQueryTests::input};
    EXPECT_EQ(48, std::distance(query.begin(), query.end()));
}

TEST(BAM_ZmwGroupQuery, whitelist_query_returns_nothing_from_empty_whitelist)
{
    const std::vector<int32_t> whitelist;

    size_t zmwCount = 0;
    size_t recordCount = 0;
    PacBio::BAM::ZmwGroupQuery query{whitelist, ZmwQueryTests::input};
    for (const auto& zmw : query) {
        ++zmwCount;
        recordCount += zmw.size();
    }
    EXPECT_EQ(0, zmwCount);
    EXPECT_EQ(0, recordCount);
}

TEST(BAM_ZmwGroupQuery, whitelist_query_returns_only_requested_zmws)
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
        recordCount += zmw.size();
    }
    EXPECT_EQ(3, zmwCount);
    EXPECT_EQ(48, recordCount);
}

TEST(BAM_ZmwGroupQuery, round_robin_query_can_return_records_applying_dataset_filter)
{
    size_t zmwCount = 0;
    size_t totalRecordCount = 0;
    std::vector<size_t> numRecordsPerZmw;
    PacBio::BAM::ZmwGroupQuery query{ZmwQueryTests::input,
                                     PacBio::BAM::ZmwFileIterationMode::SEQUENTIAL,
                                     PacBio::BAM::DataSetFilterMode::APPLY};
    for (const auto& zmw : query) {
        ++zmwCount;
        totalRecordCount += zmw.size();
        numRecordsPerZmw.push_back(zmw.size());
    }

    // zmw < 1816
    EXPECT_EQ(15, zmwCount);
    EXPECT_EQ(150, totalRecordCount);

    const std::vector<size_t> expectedNumRecordsPerZmw{2,  21, 13, 1, 5, 13, 1, 34,
                                                       12, 2,  20, 5, 3, 7,  11};
    EXPECT_TRUE(std::equal(numRecordsPerZmw.cbegin(), numRecordsPerZmw.cend(),
                           expectedNumRecordsPerZmw.cbegin()));
}

TEST(BAM_ZmwGroupQuery, round_robin_query_can_return_records_ignoring_dataset_filter)
{
    size_t zmwCount = 0;
    size_t recordCount = 0;
    std::vector<int32_t> holeNumbers;
    PacBio::BAM::ZmwGroupQuery query{ZmwQueryTests::input,
                                     PacBio::BAM::ZmwFileIterationMode::ROUND_ROBIN,
                                     PacBio::BAM::DataSetFilterMode::IGNORE};
    for (const auto& zmw : query) {
        ++zmwCount;
        if (!zmw.empty()) {
            holeNumbers.push_back(zmw.front().HoleNumber());
        }
        recordCount += zmw.size();
    }
    EXPECT_EQ(90, zmwCount);       // 30 + 30 + 30
    EXPECT_EQ(1220, recordCount);  // 432 + 409 + 379

    const std::vector<int32_t> expectedHoleNumbers{
        55,      // file 1
        54636,   // file 2
        109034,  // file 3
        480,     // file 1
        54680,   // file 2
        109043   // file 3
    };
    EXPECT_TRUE(
        std::equal(holeNumbers.cbegin(), holeNumbers.cbegin() + 6, expectedHoleNumbers.cbegin()));
}

TEST(BAM_ZmwGroupQuery, sequential_query_can_return_records_ignoring_dataset_filter)
{
    size_t zmwCount = 0;
    size_t recordCount = 0;
    std::vector<int32_t> holeNumbers;
    PacBio::BAM::ZmwGroupQuery query{ZmwQueryTests::input,
                                     PacBio::BAM::ZmwFileIterationMode::SEQUENTIAL,
                                     PacBio::BAM::DataSetFilterMode::IGNORE};
    for (const auto& zmw : query) {
        ++zmwCount;
        if (!zmw.empty() && holeNumbers.size() < 5) {
            holeNumbers.push_back(zmw.front().HoleNumber());
        }
        recordCount += zmw.size();
    }
    EXPECT_EQ(90, zmwCount);       // 30 + 30 + 30
    EXPECT_EQ(1220, recordCount);  // 432 + 409 + 379

    // first 5 ZMWs from file 1
    const std::vector<int32_t> expectedHoleNumbers{55, 480, 678, 918, 1060};
    EXPECT_TRUE(
        std::equal(holeNumbers.cbegin(), holeNumbers.cbegin() + 5, expectedHoleNumbers.cbegin()));
}

TEST(BAM_ZmwGroupQuery, can_apply_custom_pbi_filter)
{
    size_t zmwCount = 0;
    size_t recordCount = 0;
    std::vector<int32_t> holeNumbers;

    PacBio::BAM::PbiFilter zmwRange{
        {PacBio::BAM::PbiZmwFilter{1600, PacBio::BAM::Compare::GREATER_THAN_EQUAL},
         PacBio::BAM::PbiZmwFilter{1700, PacBio::BAM::Compare::LESS_THAN}}};
    PacBio::BAM::ZmwGroupQuery query{ZmwQueryTests::input, zmwRange};

    for (const auto& zmw : query) {
        ++zmwCount;
        ASSERT_TRUE(!zmw.empty());
        holeNumbers.push_back(zmw.front().HoleNumber());
        recordCount += zmw.size();
    }
    EXPECT_EQ(3, zmwCount);
    EXPECT_EQ(15, recordCount);  // 5 + 3 + 7

    const std::vector<int32_t> expectedHoleNumbers{1603, 1638, 1640};
    EXPECT_TRUE(std::equal(holeNumbers.cbegin(), holeNumbers.cend(), expectedHoleNumbers.cbegin()));
}
