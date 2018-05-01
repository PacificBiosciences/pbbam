// Author: Derek Barnett

#include <cstddef>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace FastaTests {

static void CheckSequence(const size_t index, const FastaSequence& seq)
{
    SCOPED_TRACE("checking FASTA seq:" + std::to_string(index));
    switch (index) {
        case 0:
            EXPECT_EQ("1", seq.Name());
            EXPECT_EQ(
                "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACGCAGCTCCG"
                "CCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGGCGCAGGCG",
                seq.Bases());
            break;

        case 1:
            EXPECT_EQ("2", seq.Name());
            EXPECT_EQ(
                "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACGCAGCTCCG"
                "CCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAAC",
                seq.Bases());
            break;

        case 2:
            EXPECT_EQ("3", seq.Name());
            EXPECT_EQ(
                "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACACCCTAACCCCA"
                "ACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCT",
                seq.Bases());
            break;

        default:
            ASSERT_TRUE(false);  // invalid index
    }
}

}  // namespace FastaTests

TEST(FastaSequenceTest, BasicConstructorOk)
{
    FastaSequence seq{"1", "GATTACA"};
    EXPECT_EQ("1", seq.Name());
    EXPECT_EQ("GATTACA", seq.Bases());
}

TEST(FastaReaderTest, IterableOk)
{
    const std::string fn = PbbamTestsConfig::GeneratedData_Dir + "/normal.fa";
    FastaReader reader{fn};

    size_t count = 0;
    FastaSequence seq;
    while (reader.GetNext(seq)) {
        FastaTests::CheckSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(FastaReaderTest, ReadAllOk)
{
    const std::string fn = PbbamTestsConfig::GeneratedData_Dir + "/normal.fa";

    size_t count = 0;
    for (const auto& seq : FastaReader::ReadAll(fn)) {
        FastaTests::CheckSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(FastaSequenceQueryTest, FromFastaFilename)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";

    {
        size_t count = 0;
        FastaSequenceQuery query{fn};
        for (const auto& seq : query) {
            UNUSED(seq);
            ++count;
        }
        EXPECT_EQ(1, count);
    }

    {
        FastaSequenceQuery query{fn};
        const auto first = query.cbegin();
        const auto& seq = *first;
        EXPECT_EQ("lambda_NEB3011", seq.Name());
    }
}

TEST(FastaSequenceQueryTest, FromDataSet)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/referenceset.xml";

    {
        size_t count = 0;
        FastaSequenceQuery query{fn};
        for (const auto& seq : query) {
            UNUSED(seq);
            ++count;
        }
        EXPECT_EQ(5, count);  // 1 from lambda, 4 from chimera
    }
    {
        FastaSequenceQuery query{fn};
        const auto first = query.cbegin();
        const auto& seq = *first;
        EXPECT_EQ("lambda_NEB3011", seq.Name());
    }
}
