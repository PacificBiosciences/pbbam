// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/ReadAccuracyQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(ReadAccuracyQueryTest, QueryOk)
{
    const auto bamFile = BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}};

    {
        ReadAccuracyQuery query(0.901, Compare::GREATER_THAN_EQUAL, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(4, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_GE(r.ReadAccuracy(), 0.901);
        }
        EXPECT_EQ(4, count);
    }
    {
        ReadAccuracyQuery query(0.95, Compare::GREATER_THAN_EQUAL, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(0, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_GE(r.ReadAccuracy(), 0.901);
        }
        EXPECT_EQ(0, count);
    }
}
