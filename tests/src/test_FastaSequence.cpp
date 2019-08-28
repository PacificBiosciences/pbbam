// Author: Derek Barnett

#include <cstddef>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/FastaSequence.h>

#include "FastxTests.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace FastaSequenceTests {
}  // namespace FastaSequenceTests

TEST(FastaSequenceTest, can_construct_from_seq_name_and_bases)
{
    FastaSequence seq{"1", "GATTACA"};
    EXPECT_EQ("1", seq.Name());
    EXPECT_EQ("GATTACA", seq.Bases());
}

TEST(FastaSequenceTest, can_construct_from_seq_name_and_bases_whitespaces)
{
    FastaSequence seq{"1", "GATTACA\n"};
    EXPECT_EQ("1", seq.Name());
    EXPECT_EQ("GATTACA", seq.Bases());
}