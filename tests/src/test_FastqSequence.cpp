// Author: Derek Barnett

#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>

#include "PbbamTestData.h"

#include <pbbam/BamFileMerger.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastqReader.h>
#include <pbbam/FastqSequence.h>
#include <pbbam/FastqWriter.h>

#include "FastxTests.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(FastqSequenceTest, BasicConstructorsOk)
{
    FastqSequence seq1{"1", "GATTACA", "[[[[[[["};
    EXPECT_EQ("1", seq1.Name());
    EXPECT_EQ("GATTACA", seq1.Bases());
    EXPECT_EQ("[[[[[[[", seq1.Qualities().Fastq());

    const auto quals = std::vector<uint8_t>{58, 58, 58, 58, 58, 58, 58};
    FastqSequence seq2{"1", "GATTACA", Data::QualityValues{quals}};
    EXPECT_EQ("1", seq2.Name());
    EXPECT_EQ("GATTACA", seq2.Bases());
    EXPECT_EQ("[[[[[[[", seq2.Qualities().Fastq());
}
