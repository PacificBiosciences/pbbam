// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/SubreadLengthQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(SubreadLengthQueryTest, QueryOk)
{
    const auto bamFile = BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}};

    {
        SubreadLengthQuery query(500, Compare::GREATER_THAN_EQUAL, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(3, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_GE((r.QueryEnd() - r.QueryStart()), 500);
        }
        EXPECT_EQ(3, count);
    }
    {
        SubreadLengthQuery query(1000, Compare::GREATER_THAN_EQUAL, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(2, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_GE((r.QueryEnd() - r.QueryStart()), 1000);
        }
        EXPECT_EQ(2, count);
    }
    {
        SubreadLengthQuery query(5000, Compare::GREATER_THAN_EQUAL, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(0, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_GE((r.QueryEnd() - r.QueryStart()), 5000);
        }
        EXPECT_EQ(0, count);
    }
}
