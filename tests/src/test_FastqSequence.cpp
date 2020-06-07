// Author: Derek Barnett

#include <pbbam/FastqSequence.h>

#include <gtest/gtest.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(FastqSequenceTest, BasicConstructorsOk)
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
