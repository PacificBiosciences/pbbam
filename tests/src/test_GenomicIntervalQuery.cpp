// Author: Derek Barnett

#include <iostream>
#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/GenomicIntervalQuery.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace GenomicIntervalQueryTests {
const std::string inputBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
}  // namespace GenomicIntervalQueryTests

TEST(GenomicIntervalQueryTest, ReuseQueryAndCountRecords)
{
    const std::string rname = "lambda_NEB3011";

    BamFile bamFile(GenomicIntervalQueryTests::inputBamFn);

    // setup with normal interval
    int count = 0;
    GenomicInterval interval(rname, 5000, 6000);
    GenomicIntervalQuery query(interval, bamFile);
    for (const BamRecord& record : query) {
        UNUSED(record);
        ++count;
    }
    EXPECT_EQ(2, count);

    // adjust interval and pass back in
    count = 0;
    interval.Start(9300);
    interval.Stop(9400);
    query.Interval(interval);
    for (const BamRecord& record : query) {
        UNUSED(record);
        ++count;
    }
    EXPECT_EQ(2, count);

    // adjust again (empty region)
    count = 0;
    interval.Name(rname);
    interval.Start(1000);
    interval.Stop(2000);
    query.Interval(interval);
    for (const BamRecord& record : query) {
        UNUSED(record);
        ++count;
    }
    EXPECT_EQ(0, count);

    // unknown ref
    count = 0;
    interval.Name("does not exist");
    interval.Start(0);
    interval.Stop(100);
    EXPECT_THROW(query.Interval(interval), std::runtime_error);
    for (const BamRecord& record : query) {  // iteration is still safe, just returns no data
        UNUSED(record);
        ++count;
    }
    EXPECT_EQ(0, count);

    // adjust again - make sure we can read a real region after an invalid one
    interval.Name(rname);
    interval.Start(5000);
    interval.Stop(6000);
    query.Interval(interval);
    count = 0;
    for (const BamRecord& record : query) {
        UNUSED(record);
        ++count;
    }
    EXPECT_EQ(2, count);
}

TEST(GenomicIntervalQueryTest, NonConstBamRecord)
{
    EXPECT_NO_THROW({
        BamFile bamFile(GenomicIntervalQueryTests::inputBamFn);
        int count = 0;

        GenomicInterval interval("lambda_NEB3011", 8000, 10000);
        GenomicIntervalQuery query(interval, bamFile);
        for (BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(2, count);
    });
}

TEST(GenomicIntervalQueryTest, MissingBaiShouldThrow)
{
    GenomicInterval interval("lambda_NEB3011", 0, 100);
    const std::string phi29Bam = PbbamTestsConfig::Data_Dir + "/phi29.bam";
    const std::string hasBaiBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";

    {  // single file, missing BAI
        EXPECT_THROW(GenomicIntervalQuery query(interval, phi29Bam), std::runtime_error);
    }

    {  // from dataset, all missing BAI
        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        EXPECT_THROW(GenomicIntervalQuery query(interval, ds), std::runtime_error);
    }

    {  // from dataset, mixed BAI presence
        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(
            ExternalResource("PacBio.AlignmentFile.AlignmentBamFile", hasBaiBam));
        EXPECT_THROW(GenomicIntervalQuery query(interval, ds), std::runtime_error);
    }
}
