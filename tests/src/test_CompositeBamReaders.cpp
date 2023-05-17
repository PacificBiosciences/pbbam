#include <pbbam/CompositeBamReader.h>

#include <iterator>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio::BAM;

namespace CompositeBamReaderTests {

const std::string alignedBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
const std::string aligned2BamFn = PbbamTestsConfig::Data_Dir + "/aligned2.bam";
const std::string phi29BamFn = PbbamTestsConfig::Data_Dir + "/phi29.bam";

}  // namespace CompositeBamReaderTests

TEST(BAM_GenomicIntervalCompositeBamReader, can_be_reused)
{
    const std::string refName{"lambda_NEB3011"};
    const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::alignedBamFn},
                                        BamFile{CompositeBamReaderTests::alignedBamFn}};

    // setup with normal interval
    PacBio::Data::GenomicInterval interval{refName, 5000, 6000};
    GenomicIntervalCompositeBamReader reader{interval, bamFiles};
    EXPECT_EQ(4, std::distance(reader.begin(), reader.end()));

    // adjust interval and pass back in
    interval.Start(9300);
    interval.Stop(9400);
    reader.Interval(interval);
    EXPECT_EQ(4, std::distance(reader.begin(), reader.end()));

    // adjust again (empty region)
    interval.Name(refName);
    interval.Start(1000);
    interval.Stop(2000);
    reader.Interval(interval);
    EXPECT_EQ(0, std::distance(reader.begin(), reader.end()));

    // unknown ref
    interval.Name("does not exist");
    interval.Start(0);
    interval.Stop(100);
    EXPECT_THROW(reader.Interval(interval), std::runtime_error);
    // iteration is still safe, just returns no data
    EXPECT_EQ(0, std::distance(reader.begin(), reader.end()));

    // adjust again - make sure we can read a real region after an invalid one
    interval.Name(refName);
    interval.Start(5000);
    interval.Stop(6000);
    reader.Interval(interval);
    EXPECT_EQ(4, std::distance(reader.begin(), reader.end()));
}

TEST(BAM_GenomicIntervalCompositeBamReader, throws_on_missing_bai)
{
    const PacBio::Data::GenomicInterval interval{"lambda_NEB3011", 0, 100};

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

TEST(BAM_GenomicIntervalCompositeBamReader, can_be_intialized_without_an_interval)
{
    const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::alignedBamFn},
                                        BamFile{CompositeBamReaderTests::alignedBamFn}};

    // setup without normal interval
    GenomicIntervalCompositeBamReader reader{bamFiles};
    EXPECT_EQ(0, std::distance(reader.begin(), reader.end()));

    // pass in actual interval
    PacBio::Data::GenomicInterval interval{"lambda_NEB3011", 9300, 9400};
    reader.Interval(interval);
    EXPECT_EQ(4, std::distance(reader.begin(), reader.end()));
}

// clang-format off
TEST(BAM_PbiFilterCompositeBamReader, can_handle_normal_filters)
{
    const BamFile pbiFilterBam{PbbamTestsConfig::Data_Dir + "/group/test2.bam"};

    const std::vector<BamFile> bamFiles{
        BamFile{PbbamTestsConfig::Data_Dir + "/group/test2.bam"},
        BamFile{PbbamTestsConfig::Data_Dir + "/group/test2.bam"}}; // duplicated on purpose

    {
        const std::int32_t minQLen = 500;
        PbiFilterCompositeBamReader<> reader{PbiQueryLengthFilter{minQLen, Compare::GREATER_THAN_EQUAL},
                                             bamFiles};
        EXPECT_EQ(6, reader.NumReads());
        const auto count = std::count_if(reader.begin(), reader.end(),
            [&](const BamRecord& r) {
                const auto qLen = r.QueryEnd() - r.QueryStart();
                return qLen >= minQLen;
            });
        EXPECT_EQ(6, count);
    }
    {
        // all records aligned to reverse strand && pos >= 9200
        const std::string queryName{"m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"};
        const PacBio::Data::Strand strand = PacBio::Data::Strand::REVERSE;
        const std::uint32_t minPos = 9200;
        const auto filter = PbiFilter::Intersection({
            PbiAlignedStrandFilter{strand},
            PbiReferenceStartFilter{minPos, Compare::GREATER_THAN_EQUAL}});

        PbiFilterCompositeBamReader<> reader{filter, bamFiles};
        EXPECT_EQ(2, reader.NumReads());
        const auto count = std::count_if(reader.begin(), reader.end(),
            [&](const BamRecord& r) {
                return r.AlignedStrand() == strand &&
                       r.ReferenceStart() >= static_cast<PacBio::Data::Position>(minPos) &&
                       r.FullName() == queryName;
            });
        EXPECT_EQ(2, count);
    }
    {
        // all records aligned to forward strand && pos >= 9200
        const std::string queryName{"m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2114_2531"};
        const PacBio::Data::Strand strand = PacBio::Data::Strand::FORWARD;
        const std::uint32_t minPos = 9200;
        const auto filter = PbiFilter::Intersection({
            PbiAlignedStrandFilter{strand},
            PbiReferenceStartFilter{minPos, Compare::GREATER_THAN_EQUAL}});

        PbiFilterCompositeBamReader<> reader{filter, pbiFilterBam};
        EXPECT_EQ(1, reader.NumReads());
        const auto count = std::count_if(reader.begin(), reader.end(),
            [&](const BamRecord& r) {
                return r.AlignedStrand() == strand &&
                       r.ReferenceStart() >= static_cast<PacBio::Data::Position>(minPos) &&
                       r.FullName() == queryName;
            });
        EXPECT_EQ(1, count);
    }
    {
        // all records from RG ("b89a4406") with numMatches >= 1200
        const std::string rg{"b89a4406"};
        const std::size_t minNumMatches = 1200;
        const auto filter = PbiFilter::Intersection({
            PbiReadGroupFilter{rg},
            PbiNumMatchesFilter{minNumMatches, Compare::GREATER_THAN_EQUAL}});

        PbiFilterCompositeBamReader<> reader{filter, bamFiles};
        EXPECT_EQ(4, reader.NumReads());

        std::size_t i = 0;
        const auto count = std::count_if(reader.begin(), reader.end(),
            [&](const BamRecord& r) {
                const std::string qName =
                    (i < 2 ? "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055"
                           : "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/4101_5571");
                ++i;

                return r.ReadGroupId() == rg &&
                       r.NumMatches() >= 1200 &&
                       r.FullName() == qName;
            });
        EXPECT_EQ(4, count);

    }
}
// clang-format on

TEST(BAM_SequentialCompositeBamReader, expected_record_count_across_files)
{
    const std::vector<BamFile> bamFiles{BamFile{CompositeBamReaderTests::alignedBamFn},
                                        BamFile{CompositeBamReaderTests::alignedBamFn}};

    EXPECT_NO_THROW({
        PacBio::BAM::SequentialCompositeBamReader reader{bamFiles};
        EXPECT_EQ(8, std::distance(reader.begin(), reader.end()));
    });
}
