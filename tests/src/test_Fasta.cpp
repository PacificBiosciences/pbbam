// Author: Derek Barnett

#include <cstddef>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/FastaWriter.h>
#include <pbbam/Unused.h>
#include <boost/algorithm/string.hpp>

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

TEST(FastaReaderTest, RangeForOk)
{
    const std::string fn = PbbamTestsConfig::GeneratedData_Dir + "/normal.fa";

    size_t count = 0;
    FastaReader reader{fn};
    for (const auto& seq : reader) {
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

TEST(FastaWriterTest, WriteFastaSequence)
{
    const std::string outFasta = PbbamTestsConfig::GeneratedData_Dir + "/out.fa";
    const FastaSequence seq{"name", "GATTACA"};

    {
        FastaWriter writer{outFasta};
        writer.Write(seq);
    }

    const auto seqs = FastaReader::ReadAll(outFasta);
    ASSERT_EQ(1, seqs.size());
    EXPECT_EQ(seq.Name(), seqs[0].Name());
    EXPECT_EQ(seq.Bases(), seqs[0].Bases());

    remove(outFasta.c_str());
}

TEST(FastaWriterTest, WriteBamRecord)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
    const std::string outFasta = PbbamTestsConfig::GeneratedData_Dir + "/out.fa";

    {
        FastaWriter writer{outFasta};
        EntireFileQuery query{fn};
        for (const auto& bam : query)
            writer.Write(bam);
    }
    const auto seqs = FastaReader::ReadAll(outFasta);
    ASSERT_EQ(4, seqs.size());

    const std::string name1{"singleInsertion/100/0_49"};
    const std::string name2{"singleInsertion/200/0_49"};
    const std::string name3{"singleInsertion/100/0_111"};
    const std::string name4{"singleInsertion/100/0_111"};

    EXPECT_EQ(name1, seqs[0].Name());
    EXPECT_EQ(name2, seqs[1].Name());
    EXPECT_EQ(name3, seqs[2].Name());
    EXPECT_EQ(name4, seqs[3].Name());

    const std::string bases1{"GGCTGCAGGTACAGCGGTCAGGAGGCCAATTGATGCCGGACTGGCTGAT"};
    const std::string bases2{"GGCTGCAGGTACAGCGGTCAGGAGGCCAATTGATGCCGGACTGGCTGAT"};
    const std::string bases3{
        "TTTGGCTGCAGGTACAGCGGTCAGGAGGCCAATTGATGCCGGACTGGCTGATAAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGAG"
        "CAGCACGGTAAACAGCGGCAA"};
    const std::string bases4{
        "AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGAGCAGCACGGTAAACAGCGGCAAATCAGCCAGTCCGGCATCAATTGGCCTCCTG"
        "ACCGCTGTACCTGCAGCCAAA"};

    remove(outFasta.c_str());
}

TEST(FastaWriterTest, WriteStrings)
{
    const std::string outFasta = PbbamTestsConfig::GeneratedData_Dir + "/out.fa";
    const std::string name = "name";
    const std::string bases = "GATTACA";

    {
        FastaWriter writer{outFasta};
        writer.Write(name, bases);
    }

    const auto seqs = FastaReader::ReadAll(outFasta);
    ASSERT_EQ(1, seqs.size());
    EXPECT_EQ(name, seqs[0].Name());
    EXPECT_EQ(bases, seqs[0].Bases());

    remove(outFasta.c_str());
}

TEST(FastaReaderTest, WindowsFormattedFasta)
{
    const std::string fn =
        PbbamTestsConfig::Data_Dir + "/test_windows_formatted_fasta/windows.fasta";

    {
        size_t count = 0;
        FastaReader reader{fn};
        FastaSequence seq;
        while (reader.GetNext(seq)) {
            ++count;
            bool endOK = (boost::algorithm::ends_with(seq.Name(), "5p") ||
                          boost::algorithm::ends_with(seq.Name(), "3p"));
            EXPECT_TRUE(endOK);
        }
        EXPECT_EQ(7, count);  // 7 primers in total
    }
}
