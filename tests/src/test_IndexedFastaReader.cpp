// Author: Derek Barnett

#include <iostream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/IndexedFastaReader.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace IndexedFastaReaderTests {

const std::string lambdaFasta = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";
const std::string singleInsertionBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";

}  // namespace IndexedFastaReaderTests

TEST(IndexedFastaReaderTests, PrintSingleInsertion)
{
    IndexedFastaReader r(IndexedFastaReaderTests::lambdaFasta);

    // Open BAM file
    BamFile bamFile(IndexedFastaReaderTests::singleInsertionBam);
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

    // {
    //     std::ostringstream output;
    //     auto itSS = bamQuery.begin();
    //     {
    //         const auto recordSS = *itSS;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true, true) << std::endl;
    //         output << std::endl;
    //     }
    //     ++itSS;
    //     {
    //         const auto recordSS = *itSS;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true, true) << std::endl;
    //         output << std::endl;
    //     }
    //     ++itSS;
    //     {
    //         const auto recordSS = *itSS;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true, true) << std::endl;
    //         output << std::endl;
    //     }
    //     ++itSS;
    //     {
    //         const auto recordSS = *itSS;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::GENOMIC, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::GENOMIC, true, true) << std::endl;
    //         output << std::endl;
    //         output << r.ReferenceSubsequence(recordSS, Orientation::NATIVE, true, true) << std::endl;
    //         output << recordSS.Sequence(Orientation::NATIVE, true, true) << std::endl;
    //     }
    //     std::cerr << output.str();
    // }
}

TEST(IndexedFastaReaderTests, ReadLambda)
{
    IndexedFastaReader r(IndexedFastaReaderTests::lambdaFasta);

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

    // bad subsequence
}

TEST(IndexedFastaReaderTests, Errors)
{
    IndexedFastaReader r(IndexedFastaReaderTests::lambdaFasta);

    //
    // attempt access without "opening"
    //
    // EXPECT_THROW(r.NumSequences(), std::exception);
    // EXPECT_THROW(r.HasSequence("lambda_NEB3011"), std::exception);
    // EXPECT_THROW(r.SequenceLength("lambda_NEB3011"), std::exception);
    // EXPECT_THROW(r.Subsequence("lambda_NEB3011:0-10"), std::exception);

    //
    // invalid accesses after opening
    //
    EXPECT_THROW(r.SequenceLength("dog"), std::exception);
    EXPECT_THROW(r.Subsequence("dog:0-10"), std::exception);
}

TEST(IndexedFastaReaderTests, Names)
{
    IndexedFastaReader r(IndexedFastaReaderTests::lambdaFasta);
    std::vector<std::string> names = {"lambda_NEB3011"};

    // Test all-name request
    EXPECT_EQ(names, r.Names());

    // Test single-name query
    EXPECT_EQ(names[0], r.Name(0));

    // invalid name acces (out of range)
    EXPECT_THROW(r.Name(1), std::exception);
}
