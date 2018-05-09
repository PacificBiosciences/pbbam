// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/ZmwQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;

//TEST(EntireFileQueryTest, CountRecords)
//{
//    EXPECT_NO_THROW(
//    {
//        // open input BAM file
//        BamFile bamFile(inputBamFn);

//        // count records
//        int count = 0;
//        EntireFileQuery entireFile(bamFile);
//        for (const BamRecord& record : entireFile) {
//            ()record;
//            ++count;
//        }

//        EXPECT_EQ(3307, count);
//    });
//}
