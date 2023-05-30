#include <pbbam/FastaWriter.h>

#include <cstddef>

#include <gtest/gtest.h>

#include <boost/algorithm/string.hpp>

#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastaSequenceQuery.h>

#include "FastxTests.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace FastaWriterTests {}  // namespace FastaWriterTests

TEST(BAM_FastaWriter, throws_on_empty_filename)
{
    EXPECT_THROW(FastaWriter writer{""}, std::runtime_error);
}

TEST(BAM_FastaWriter, throws_on_invalid_extension)
{
    EXPECT_THROW(FastaWriter writer{"wrong.ext"}, std::runtime_error);
}

TEST(BAM_FastaWriter, can_write_fasta_sequence)
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

    std::remove(outFasta.c_str());
}

TEST(BAM_FastaWriter, can_write_fasta_from_bam)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
    const std::string outFasta = PbbamTestsConfig::GeneratedData_Dir + "/out.fa";

    {
        FastaWriter writer{outFasta};
        EntireFileQuery query{fn};
        for (const auto& bam : query) {
            writer.Write(bam);
        }
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

    std::remove(outFasta.c_str());
}

TEST(BAM_FastaWriter, can_write_fasta_from_strings)
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

    std::remove(outFasta.c_str());
}
