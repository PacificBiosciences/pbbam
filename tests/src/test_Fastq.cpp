// Author: Derek Barnett

#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>

#include "PbbamTestData.h"

#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastqReader.h>
#include <pbbam/FastqSequence.h>
#include <pbbam/FastqWriter.h>

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

TEST(FastqWriterTest, WriteFastqSequence)
{
    const std::string outFastq = PbbamTestsConfig::GeneratedData_Dir + "/out.fq";
    const FastqSequence seq{"name", "GATTACA", "!!!!!!!"};

    {
        FastqWriter writer{outFastq};
        writer.Write(seq);
    }

    const auto seqs = FastqReader::ReadAll(outFastq);
    ASSERT_EQ(1, seqs.size());
    EXPECT_EQ(seq.Name(), seqs[0].Name());
    EXPECT_EQ(seq.Bases(), seqs[0].Bases());
    EXPECT_EQ(seq.Qualities(), seqs[0].Qualities());

    remove(outFastq.c_str());
}

TEST(FastqWriterTest, WriteBamRecord)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/unmap1.bam";
    const std::string outFastq = PbbamTestsConfig::GeneratedData_Dir + "/out.fq";

    {
        FastqWriter writer{outFastq};
        EntireFileQuery query{fn};
        for (const auto& bam : query)
            writer.Write(bam);
    }
    const auto seqs = FastqReader::ReadAll(outFastq);
    ASSERT_EQ(1, seqs.size());

    const std::string name{"test/1/0_100"};
    const std::string bases{
        "GATCGCACTGAAAATCTGGATATAGAACGTGTGCAAATGATTGTCTCTACCGTTCCGTAAAAATTATTGCTAATTAGCAATGATTTTAAG"
        "CTAATTAGTT"};
    const std::string quals{
        "CCCCCCCCCCCCCCCCCCCACCCCCACCCCCCCCCCCCB;CCCAACCCCCCCCCCCCCD=B9BCABCBCB>BBBC@B<<@BA;BCC?B>"
        "A<<@(?:4==4"};

    EXPECT_EQ(name, seqs[0].Name());
    EXPECT_EQ(bases, seqs[0].Bases());
    EXPECT_EQ(quals, seqs[0].Qualities().Fastq());

    remove(outFastq.c_str());
}

TEST(FastqWriterTest, WriteStrings)
{
    const std::string outFastq = PbbamTestsConfig::GeneratedData_Dir + "/out.fq";
    const std::string name = "name";
    const std::string bases = "GATTACA";
    const std::string quals = "!!!!!!!";

    {
        FastqWriter writer{outFastq};
        writer.Write(name, bases, quals);
    }

    const auto seqs = FastqReader::ReadAll(outFastq);
    ASSERT_EQ(1, seqs.size());
    EXPECT_EQ(name, seqs[0].Name());
    EXPECT_EQ(bases, seqs[0].Bases());

    remove(outFastq.c_str());
}
