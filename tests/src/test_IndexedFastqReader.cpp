// Author: Derek Barnett

#include <pbbam/IndexedFastqReader.h>

#include <gtest/gtest.h>

#include <htslib/hts.h>

#include "FastxTests.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(IndexedFastqReaderTest, throws_on_empty_filename)
{
    EXPECT_THROW(IndexedFastqReader reader{""}, std::runtime_error);
}

TEST(IndexedFastqReaderTest, throws_on_invalid_extension)
{
    EXPECT_THROW(IndexedFastqReader reader{"wrong.ext"}, std::runtime_error);
}

TEST(IndexedFastqReaderTest, can_open_text_fastq)
{
    const auto& fn = FastxTests::simpleFastqFn;
    EXPECT_NO_THROW(IndexedFastqReader reader{fn});
}

TEST(IndexedFastqReaderTest, throws_on_gzip_fastq)
{
    const auto& fn = FastxTests::simpleFastqGzipFn;
    EXPECT_THROW(IndexedFastqReader reader{fn}, std::runtime_error);
}

TEST(IndexedFastqReaderTest, can_open_bgzf_fastq_for_reading)
{
    const auto& fn = FastxTests::simpleFastqBgzfFn;
    EXPECT_NO_THROW(IndexedFastqReader reader{fn});
}

TEST(IndexedFastqReaderTest, can_query_index_for_metadata)
{
    const IndexedFastqReader r{FastxTests::simpleFastqFn};

    EXPECT_TRUE(r.HasSequence("seq1"));
    EXPECT_FALSE(r.HasSequence("nope"));

    EXPECT_EQ(8, r.NumSequences());
    EXPECT_EQ(63, r.SequenceLength("seq5"));

    const auto& names = r.Names();
    ASSERT_EQ(FastxTests::ExpectedFastq.size(), names.size());
    for (size_t i = 0; i < names.size(); ++i)
        EXPECT_EQ(FastxTests::ExpectedFastq.at(i).Name(), names.at(i));
}

TEST(IndexedFastqReaderTest, subsequence_from_text_fastq)
{
    IndexedFastqReader r{FastxTests::simpleFastqFn};
    {
        const std::string expectedSeq{"GCATGCATGC"};
        const QualityValues expectedQuals{"~}|{zyxwvu"};

        const auto subsequence = r.Subsequence("seq2", 0, 10);

        EXPECT_EQ(expectedSeq, subsequence.first);
        EXPECT_EQ(expectedQuals, subsequence.second.Fastq());
    }
    {
        const std::string expectedSeq{"ATGCATGCAT"};
        const QualityValues expectedQuals{R"(`_^]\[ZYXW)"};

        const auto subsequence = r.Subsequence("seq6", 30, 40);

        EXPECT_EQ(expectedSeq, subsequence.first);
        EXPECT_EQ(expectedQuals, subsequence.second.Fastq());
    }
}

TEST(IndexedFastqReaderTest, subsequence_from_bgzf_fastq)
{
    IndexedFastqReader r{FastxTests::simpleFastqBgzfFn};
    {
        const std::string expectedSeq{"GCATGCATGC"};
        const QualityValues expectedQuals{"~}|{zyxwvu"};

        const auto subsequence = r.Subsequence("seq2", 0, 10);

        EXPECT_EQ(expectedSeq, subsequence.first);
        EXPECT_EQ(expectedQuals, subsequence.second.Fastq());
    }
    {
        const std::string expectedSeq{"ATGCATGCAT"};
        const QualityValues expectedQuals{R"(`_^]\[ZYXW)"};

        const auto subsequence = r.Subsequence("seq6", 30, 40);

        EXPECT_EQ(expectedSeq, subsequence.first);
        EXPECT_EQ(expectedQuals, subsequence.second.Fastq());
    }
}

TEST(IndexedFastqReaderTest, returns_empty_result_from_empty_region)
{
    IndexedFastqReader r{FastxTests::simpleFastqFn};
    const auto subsequence = r.Subsequence("seq2", 0, 0);
    EXPECT_TRUE(subsequence.first.empty());
    EXPECT_TRUE(subsequence.second.empty());
}

TEST(IndexedFastqReaderTest, throws_if_region_is_malformated)
{
    IndexedFastqReader r{FastxTests::simpleFastqFn};

    // start > end
    EXPECT_THROW(r.Subsequence("seq2", 10, 5), std::runtime_error);

    // start, end < 0
    EXPECT_THROW(r.Subsequence("seq2", -1, 5), std::runtime_error);
    EXPECT_THROW(r.Subsequence("seq2", 5, -1), std::runtime_error);
    EXPECT_THROW(r.Subsequence("seq2", -2, -1), std::runtime_error);
}

TEST(IndexedFastqReaderTest, returns_available_length_if_region_is_longer)
{
    // i.e. like std::string::substr()

    IndexedFastqReader r{FastxTests::simpleFastqFn};
    const auto subsequence = r.Subsequence("seq2", 0, 1000);
    EXPECT_EQ(63, subsequence.first.size());
}
