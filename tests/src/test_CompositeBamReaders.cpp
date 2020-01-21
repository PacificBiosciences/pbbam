// Author: Derek Barnett

#include <tuple>

#include <gtest/gtest.h>

#include <pbbam/CompositeBamReader.h>

#include "PbbamTestData.h"

using namespace PacBio::BAM;

namespace CompositeBamReaderTests {

const std::string alignedBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
const std::string aligned2BamFn = PbbamTestsConfig::Data_Dir + "/aligned2.bam";
const std::string phi29BamFn = PbbamTestsConfig::Data_Dir + "/phi29.bam";

}  // namespace CompositeBamReaderTests

TEST(GenomicIntervalCompositeBamReaderTest, ReuseReader)
{
    const std::string refName{"lambda_NEB3011"};
    const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::alignedBamFn},
                                        BamFile{CompositeBamReaderTests::alignedBamFn}};

    // setup with normal interval
    int count = 0;
    GenomicInterval interval(refName, 5000, 6000);
    GenomicIntervalCompositeBamReader reader{interval, bamFiles};
    for (const auto& record : reader) {
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(4, count);

    // adjust interval and pass back in
    count = 0;
    interval.Start(9300);
    interval.Stop(9400);
    reader.Interval(interval);
    for (const auto& record : reader) {
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(4, count);

    // adjust again (empty region)
    count = 0;
    interval.Name(refName);
    interval.Start(1000);
    interval.Stop(2000);
    reader.Interval(interval);
    for (const auto& record : reader) {
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(0, count);

    // unknown ref
    count = 0;
    interval.Name("does not exist");
    interval.Start(0);
    interval.Stop(100);
    EXPECT_THROW(reader.Interval(interval), std::runtime_error);
    for (const auto& record : reader) {  // iteration is still safe, just returns no data
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(0, count);

    // adjust again - make sure we can read a real region after an invalid one
    interval.Name(refName);
    interval.Start(5000);
    interval.Stop(6000);
    reader.Interval(interval);
    count = 0;
    for (const auto& record : reader) {
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(4, count);
}

TEST(GenomicIntervalCompositeBamReaderTest, MissingBaiShouldThrow)
{
    const GenomicInterval interval{"lambda_NEB3011", 0, 100};

    {  // single file, missing BAI
        const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::phi29BamFn}};
        EXPECT_THROW(GenomicIntervalCompositeBamReader reader(interval, bamFiles),
                     std::runtime_error);
    }

    {  // from dataset, all missing BAI
        const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::phi29BamFn},
                                            BamFile{CompositeBamReaderTests::phi29BamFn}};
        EXPECT_THROW(GenomicIntervalCompositeBamReader reader(interval, bamFiles),
                     std::runtime_error);
    }

    {  // from dataset, mixed BAI presence
        const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::phi29BamFn},
                                            BamFile{CompositeBamReaderTests::alignedBamFn}};
        EXPECT_THROW(GenomicIntervalCompositeBamReader reader(interval, bamFiles),
                     std::runtime_error);
    }
}

TEST(GenomicIntervalCompositeBamReaderTest, InitializeWithoutInterval)
{
    const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::alignedBamFn},
                                        BamFile{CompositeBamReaderTests::alignedBamFn}};

    // setup without normal interval
    int count = 0;
    GenomicIntervalCompositeBamReader reader{bamFiles};
    for (const auto& record : reader) {
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(0, count);

    // pass in actual interval
    count = 0;
    GenomicInterval interval{"lambda_NEB3011", 9300, 9400};
    reader.Interval(interval);
    for (const auto& record : reader) {
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(4, count);
}

TEST(PbiFilterCompositeBamReaderTest, BasicFilteringOk)
{
    const auto pbiFilterBam = BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}};

    const std::vector<BamFile> bamFiles{
        BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}},
        BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}}};

    {
        PbiFilterCompositeBamReader<> reader{PbiQueryLengthFilter{500, Compare::GREATER_THAN_EQUAL},
                                             bamFiles};
        const auto numReads = reader.NumReads();
        EXPECT_EQ(6, numReads);

        int count = 0;
        for (const auto& r : reader) {
            ++count;
            EXPECT_GE((r.QueryEnd() - r.QueryStart()), 500);
        }
        EXPECT_EQ(6, count);
    }
    {
        // all records aligned to reverse strand && pos >= 9200
        const auto filter =
            PbiFilter::Intersection({PbiAlignedStrandFilter{Strand::REVERSE},
                                     PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}});

        PbiFilterCompositeBamReader<> reader{filter, bamFiles};
        const auto numReads = reader.NumReads();
        EXPECT_EQ(2, numReads);

        int count = 0;
        for (const auto& r : reader) {
            ++count;
            EXPECT_EQ(Strand::REVERSE, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(
                std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/"
                            "5615_6237"),
                r.FullName());
        }
        EXPECT_EQ(2, count);
    }
    {
        // all records aligned to forward strand && pos >= 9200
        const auto filter =
            PbiFilter::Intersection({PbiAlignedStrandFilter{Strand::FORWARD},
                                     PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}});

        PbiFilterCompositeBamReader<> reader{filter, pbiFilterBam};
        const auto numReads = reader.NumReads();
        EXPECT_EQ(1, numReads);

        int count = 0;
        for (const auto& r : reader) {
            ++count;
            EXPECT_EQ(Strand::FORWARD, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(
                std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/"
                            "2114_2531"),
                r.FullName());
        }
        EXPECT_EQ(1, count);
    }
    {
        // all records from RG ("b89a4406") with numMatches >= 1200
        const auto filter =
            PbiFilter::Intersection({PbiReadGroupFilter{"b89a4406"},
                                     PbiNumMatchesFilter{1200, Compare::GREATER_THAN_EQUAL}});

        PbiFilterCompositeBamReader<> reader{filter, bamFiles};
        const auto numReads = reader.NumReads();
        EXPECT_EQ(4, numReads);

        int count = 0;
        for (const auto& r : reader) {
            ++count;
            EXPECT_EQ(std::string("b89a4406"), r.ReadGroupId());
            EXPECT_GE((r.NumMatches()), 1200);
            if (count == 1 || count == 2)  // file1, file2
                EXPECT_EQ(
                    std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/"
                                "14743/2579_4055"),
                    r.FullName());
            else {
                if (count == 3 || count == 4) {  // file 1, file 2
                    EXPECT_EQ(
                        std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_"
                                    "X0/14743/4101_5571"),
                        r.FullName());
                }
            }
        }
        EXPECT_EQ(4, count);
    }
}

TEST(SequentialCompositeBamReaderTest, CountRecords)
{
    const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::alignedBamFn},
                                        BamFile{CompositeBamReaderTests::alignedBamFn}};

    EXPECT_NO_THROW({
        int count = 0;
        PacBio::BAM::SequentialCompositeBamReader reader{bamFiles};
        for (const auto& record : reader) {
            std::ignore = record;
            ++count;
        }
        EXPECT_EQ(8, count);
    });
}
