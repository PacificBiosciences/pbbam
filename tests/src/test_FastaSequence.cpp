// Author: Derek Barnett

#include <pbbam/FastaSequence.h>

#include <gtest/gtest.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(FastaSequenceTest, can_construct_from_seq_name_and_bases)
{
    const FastaSequence seq{"1", "GATTACA"};
    EXPECT_EQ("1", seq.Name());
    EXPECT_EQ("GATTACA", seq.Bases());
}

TEST(FastaSequenceTest, can_construct_from_seq_name_and_bases_whitespaces)
{
    const FastaSequence seq{"1", "GATTACA\n"};
    EXPECT_EQ("1", seq.Name());
    EXPECT_EQ("GATTACA", seq.Bases());
}
