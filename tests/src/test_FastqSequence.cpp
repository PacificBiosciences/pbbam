// Author: Derek Barnett

#include <pbbam/FastqSequence.h>

#include <gtest/gtest.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_FastqSequence, BasicConstructorsOk)
{
    const FastqSequence seq1{"1", "GATTACA", "[[[[[[["};
    EXPECT_EQ("1", seq1.Name());
    EXPECT_EQ("GATTACA", seq1.Bases());
    EXPECT_EQ("[[[[[[[", seq1.Qualities().Fastq());

    const std::vector<uint8_t> quals{58, 58, 58, 58, 58, 58, 58};
    const FastqSequence seq2{"1", "GATTACA", QualityValues{quals}};
    EXPECT_EQ("1", seq2.Name());
    EXPECT_EQ("GATTACA", seq2.Bases());
    EXPECT_EQ("[[[[[[[", seq2.Qualities().Fastq());
}

TEST(BAM_FastqSequence, average_base_quality_throws_on_empty_sequence)
{
    const FastqSequence fastq;
    EXPECT_THROW(fastq.AverageBaseQuality(), std::runtime_error);
}

TEST(BAM_FastqSequence, can_calculate_average_base_quality)
{
    {
        const std::vector<uint8_t> qvs{20};
        const FastqSequence fastq{"seq1", "G", Data::QualityValues{qvs}};
        EXPECT_FLOAT_EQ(fastq.AverageBaseQuality(), 20.0f);
    }
    {
        const std::vector<uint8_t> qvs{20, 20, 30, 30, 20, 20, 30, 30};
        const FastqSequence fastq("seq2", "GATTACA", Data::QualityValues{qvs});
        EXPECT_FLOAT_EQ(fastq.AverageBaseQuality(), 25.0f);
    }
    {
        const std::vector<uint8_t> qvs{40, 40, 40, 40, 40, 40, 40, 40};
        const FastqSequence fastq{"seq3", "GATTACA", Data::QualityValues{qvs}};
        EXPECT_FLOAT_EQ(fastq.AverageBaseQuality(), 40.0f);
    }
}
