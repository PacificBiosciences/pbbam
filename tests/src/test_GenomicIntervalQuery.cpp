// Author: Derek Barnett

#include <algorithm>
#include <functional>
#include <iterator>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/BaiIndexCache.h>
#include <pbbam/GenomicIntervalQuery.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace GenomicIntervalQueryTests {
const std::string inputBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
const std::string inputBamFn_2 = PbbamTestsConfig::Data_Dir + "/aligned2.bam";
}  // namespace GenomicIntervalQueryTests

TEST(GenomicIntervalQueryTest, ReuseQueryAndCountRecords)
{
    const std::string rname{"lambda_NEB3011"};
    const BamFile bamFile{GenomicIntervalQueryTests::inputBamFn};

    // setup with normal interval
    GenomicInterval interval{rname, 5000, 6000};
    GenomicIntervalQuery query{interval, bamFile};
    EXPECT_EQ(2, std::distance(query.begin(), query.end()));

    // adjust interval and pass back in
    interval.Start(9300);
    interval.Stop(9400);
    query.Interval(interval);
    EXPECT_EQ(2, std::distance(query.begin(), query.end()));

    // adjust again (empty region)
    interval.Name(rname);
    interval.Start(1000);
    interval.Stop(2000);
    query.Interval(interval);
    EXPECT_EQ(0, std::distance(query.begin(), query.end()));

    // unknown ref
    interval.Name("does not exist");
    interval.Start(0);
    interval.Stop(100);
    EXPECT_THROW(query.Interval(interval), std::runtime_error);
    // iteration is still safe, just returns no data
    EXPECT_EQ(0, std::distance(query.begin(), query.end()));

    // adjust again - make sure we can read a real region after an invalid one
    interval.Name(rname);
    interval.Start(5000);
    interval.Stop(6000);
    query.Interval(interval);
    EXPECT_EQ(2, std::distance(query.begin(), query.end()));
}

TEST(GenomicIntervalQueryTest, NonConstBamRecord)
{
    const BamFile bamFile{GenomicIntervalQueryTests::inputBamFn};
    const GenomicInterval interval{"lambda_NEB3011", 8000, 10000};
    GenomicIntervalQuery query{interval, bamFile};
    EXPECT_EQ(2, std::distance(query.begin(), query.end()));
}

TEST(GenomicIntervalQueryTest, MissingBaiShouldThrow)
{
    const GenomicInterval interval{"lambda_NEB3011", 0, 100};
    const std::string phi29Bam = PbbamTestsConfig::Data_Dir + "/phi29.bam";
    const std::string hasBaiBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";

    {  // single file, missing BAI
        EXPECT_THROW(GenomicIntervalQuery query(interval, phi29Bam), std::runtime_error);
    }

    {  // from dataset, all missing BAI
        DataSet ds;
        ds.ExternalResources().Add(ExternalResource{"PacBio.SubreadFile.SubreadBamFile", phi29Bam});
        ds.ExternalResources().Add(ExternalResource{"PacBio.SubreadFile.SubreadBamFile", phi29Bam});
        EXPECT_THROW(GenomicIntervalQuery query(interval, ds), std::runtime_error);
    }

    {  // from dataset, mixed BAI presence
        DataSet ds;
        ds.ExternalResources().Add(ExternalResource{"PacBio.SubreadFile.SubreadBamFile", phi29Bam});
        ds.ExternalResources().Add(
            ExternalResource{"PacBio.AlignmentFile.AlignmentBamFile", hasBaiBam});
        EXPECT_THROW(GenomicIntervalQuery query(interval, ds), std::runtime_error);
    }
}

TEST(GenomicIntervalQueryTest, InitializeWithoutInterval)
{
    const std::string rname = "lambda_NEB3011";

    const BamFile bamFile{GenomicIntervalQueryTests::inputBamFn};

    // setup without normal interval
    GenomicIntervalQuery query{bamFile};
    EXPECT_EQ(0, std::distance(query.begin(), query.end()));

    // pass in actual interval
    GenomicInterval interval{"lambda_NEB3011", 9300, 9400};
    query.Interval(interval);
    EXPECT_EQ(2, std::distance(query.begin(), query.end()));
}

TEST(GenomicIntervalQueryTest, CanReuseBaiIndexCache)
{
    const std::string refName{"lambda_NEB3011"};
    const std::vector<std::string> filenames{GenomicIntervalQueryTests::inputBamFn,
                                             GenomicIntervalQueryTests::inputBamFn_2};

    const DataSet ds{filenames};
    const auto indexCache = MakeBaiIndexCache(ds);

    auto checkInterval = [](GenomicIntervalQuery& query, const GenomicInterval& interval,
                            const size_t expectedCount) {
        // update query
        query.Interval(interval);

        // checkout results
        std::vector<Position> startPositions;
        for (const BamRecord& r : query) {
            EXPECT_EQ(interval.Name(), r.ReferenceName());
            EXPECT_TRUE(r.ReferenceStart() < interval.Stop());
            EXPECT_TRUE(r.ReferenceEnd() >= interval.Start());
            startPositions.push_back(r.ReferenceStart());
        }
        EXPECT_EQ(expectedCount, startPositions.size());
        EXPECT_TRUE(std::is_sorted(startPositions.cbegin(), startPositions.cend()));
    };

    // reuse cache between interval updates
    GenomicIntervalQuery query{ds, indexCache};
    {
        const GenomicInterval interval{refName, 5000, 8000};
        const size_t expectedCount = 7;
        checkInterval(query, interval, expectedCount);
    }
    {
        const GenomicInterval interval{refName, 0, 100};
        const size_t expectedCount = 1;
        checkInterval(query, interval, expectedCount);
    }
    {
        const GenomicInterval interval{refName, 9300, 9400};
        const size_t expectedCount = 2;
        checkInterval(query, interval, expectedCount);
    }

    // reuse cache in independent query
    GenomicIntervalQuery query2{ds, indexCache};
    const GenomicInterval interval{refName, 5000, 8000};
    const size_t expectedCount = 7;
    checkInterval(query2, interval, expectedCount);
}
