// Author: Derek Barnett

#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>

#include "PbbamTestData.h"

#include <pbbam/FastqReader.h>
#include <pbbam/FastqSequence.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace FastqTests {

static void CheckSequence(const size_t index, const FastqSequence& seq)
{
    SCOPED_TRACE("checking Fastq seq:" + std::to_string(index));
    switch (index) {
        case 0:
            EXPECT_EQ("1", seq.Name());
            EXPECT_EQ(
                "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
                "ACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGA"
                "GGAGAACGCAACTCCGCCGGCGCAGGCG",
                seq.Bases());
            EXPECT_EQ(
                "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                "[[[[[[[[[[[[[[[[[[[[[[[[[[[[",
                seq.Qualities().Fastq());
            break;

        case 1:
            EXPECT_EQ("2", seq.Name());
            EXPECT_EQ(
                "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
                "ACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGA"
                "GGAGAACGCAAC",
                seq.Bases());
            EXPECT_EQ(
                "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                "[[[[[[[[[[[[",
                seq.Qualities().Fastq());
            break;

        case 2:
            EXPECT_EQ("3", seq.Name());
            EXPECT_EQ(
                "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
                "ACCCTAACCCTAACACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA"
                "ACCCTAACCCCTAACCCTAACCCT",
                seq.Bases());
            EXPECT_EQ(
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                "]]]]]]]]]]]]]]]]]]]]]]]]",
                seq.Qualities().Fastq());
            break;

        default:
            ASSERT_TRUE(false);  // invalid index
    }
}

}  // namespace FastqTests

TEST(FastqSequenceTest, BasicConstructorsOk)
{
    FastqSequence seq1{"1", "GATTACA", "[[[[[[["};
    EXPECT_EQ("1", seq1.Name());
    EXPECT_EQ("GATTACA", seq1.Bases());
    EXPECT_EQ("[[[[[[[", seq1.Qualities().Fastq());

    const auto quals = std::vector<uint8_t>{58, 58, 58, 58, 58, 58, 58};
    FastqSequence seq2{"1", "GATTACA", QualityValues{quals}};
    EXPECT_EQ("1", seq2.Name());
    EXPECT_EQ("GATTACA", seq2.Bases());
    EXPECT_EQ("[[[[[[[", seq2.Qualities().Fastq());
}

TEST(FastqReaderTest, IterableOk)
{
    const std::string fn = PbbamTestsConfig::GeneratedData_Dir + "/normal.fq";
    FastqReader reader{fn};

    size_t count = 0;
    FastqSequence seq;
    while (reader.GetNext(seq)) {
        FastqTests::CheckSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(FastqReaderTest, ReadAllOk)
{
    const std::string fn = PbbamTestsConfig::GeneratedData_Dir + "/normal.fq";

    size_t count = 0;
    for (const auto& seq : FastqReader::ReadAll(fn)) {
        FastqTests::CheckSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(3, count);
}
