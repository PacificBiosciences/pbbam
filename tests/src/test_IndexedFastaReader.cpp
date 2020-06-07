// Author: Derek Barnett

#include <pbbam/IndexedFastaReader.h>

#include <string>

#include <gtest/gtest.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/EntireFileQuery.h>

#include "FastxTests.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace IndexedFastaReaderTests {

const std::string lambdaFasta = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";
const std::string singleInsertionBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";

}  // namespace IndexedFastaReaderTests

TEST(IndexedFastaReaderTest, throws_on_empty_filename)
{
    EXPECT_THROW(IndexedFastaReader reader{""}, std::runtime_error);
}

TEST(IndexedFastaReaderTest, throws_on_invalid_extension)
{
    EXPECT_THROW(IndexedFastaReader reader{"wrong.ext"}, std::runtime_error);
}

TEST(IndexedFastaReaderTest, can_open_text_fasta_for_reading)
{
    const auto& fn = FastxTests::simpleFastaFn;
    EXPECT_NO_THROW(IndexedFastaReader reader{fn});
}

TEST(IndexedFastaReaderTest, throws_on_gzip_fasta)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    EXPECT_THROW(IndexedFastaReader reader{fn}, std::runtime_error);
}

TEST(IndexedFastaReaderTest, can_open_bgzf_fasta_for_reading)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    EXPECT_NO_THROW(IndexedFastaReader reader{fn});
}

TEST(IndexedFastaReaderTest, can_fetch_subsequence_from_lambda)
{
    IndexedFastaReader r{IndexedFastaReaderTests::lambdaFasta};

    EXPECT_TRUE(r.HasSequence("lambda_NEB3011"));
    EXPECT_FALSE(r.HasSequence("dog"));
    EXPECT_EQ(1, r.NumSequences());
    EXPECT_EQ(48502, r.SequenceLength("lambda_NEB3011"));

    std::string seq = r.Subsequence("lambda_NEB3011:0-10");
    EXPECT_EQ("GGGCGGCGAC", seq);

    std::string seq2 = r.Subsequence("lambda_NEB3011", 0, 10);
    EXPECT_EQ("GGGCGGCGAC", seq2);

    // subsequence extending beyond bounds returns clipped
    std::string seq3 = r.Subsequence("lambda_NEB3011", 48400, 48600);
    EXPECT_EQ(102, seq3.length());

    // empty subsequence
    std::string emptySeq = r.Subsequence("lambda_NEB3011", 10, 10);
    EXPECT_EQ("", emptySeq);
}

TEST(IndexedFastaReaderTest, prints_clipped_and_gapped_subsequences_from_lambda)
{
    IndexedFastaReader r{IndexedFastaReaderTests::lambdaFasta};

    // Open BAM file
    const BamFile bamFile{IndexedFastaReaderTests::singleInsertionBam};
    EntireFileQuery bamQuery(bamFile);

    auto it = bamQuery.begin();
    auto record = *it++;
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::NATIVE, true));
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::NATIVE, true, true));
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::GENOMIC, true));
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::GENOMIC, true, true));
    record = *it++;
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::NATIVE, true));
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::NATIVE, true, true));
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::GENOMIC, true));
    EXPECT_EQ("GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT",
              r.ReferenceSubsequence(record, Orientation::GENOMIC, true, true));
    record = *it++;
    EXPECT_EQ(
        "----------------------------------------------------"
        "AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
        r.ReferenceSubsequence(record, Orientation::NATIVE, true));
    EXPECT_EQ("AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
              r.ReferenceSubsequence(record, Orientation::NATIVE, true, true));
    EXPECT_EQ(
        "----------------------------------------------------"
        "AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
        r.ReferenceSubsequence(record, Orientation::GENOMIC, true));
    EXPECT_EQ("AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
              r.ReferenceSubsequence(record, Orientation::GENOMIC, true, true));
    record = *it++;
    EXPECT_EQ(
        "AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA-----------------------------"
        "-----------------------",
        r.ReferenceSubsequence(record, Orientation::GENOMIC, true));
    EXPECT_EQ(
        "----------------------------------------------------TTGCCGCTGTT-"
        "ACCGTGCTGCGATCTTCTGCCATCGACGGACGTCCCACATTGGTGACTT",
        r.ReferenceSubsequence(record, Orientation::NATIVE, true));
    EXPECT_EQ("AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
              r.ReferenceSubsequence(record, Orientation::GENOMIC, true, true));
    EXPECT_EQ("TTGCCGCTGTT-ACCGTGCTGCGATCTTCTGCCATCGACGGACGTCCCACATTGGTGACTT",
              r.ReferenceSubsequence(record, Orientation::NATIVE, true, true));
}

// Come back
TEST(IndexedFastaReaderTest, throws_on_invalid_subsequence_requests)
{
    IndexedFastaReader r{IndexedFastaReaderTests::lambdaFasta};
    EXPECT_THROW(r.SequenceLength("dog"), std::exception);
    EXPECT_THROW(r.Subsequence("dog:0-10"), std::exception);
}

//
TEST(IndexedFastaReaderTest, can_fetch_name_info_from_lambda)
{
    IndexedFastaReader r{IndexedFastaReaderTests::lambdaFasta};
    const std::vector<std::string> names{"lambda_NEB3011"};

    // Test all-name request
    EXPECT_EQ(names, r.Names());

    // Test single-name query
    EXPECT_EQ(names[0], r.Name(0));

    // invalid name acces (out of range)
    EXPECT_THROW(r.Name(1), std::exception);
}
