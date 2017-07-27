// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#include <iostream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#ifdef PBBAM_TESTING
#define private public
#endif

#include "TestData.h"

#include "pbbam/IndexedFastaReader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/BamFile.h"
#include "pbbam/EntireFileQuery.h"

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string lambdaFasta = tests::Data_Dir + "/lambdaNEB.fa";
const string singleInsertionBam = tests::Data_Dir + "/aligned.bam";

TEST(IndexedFastaReaderTests, PrintSingleInsertion)
{
    IndexedFastaReader r(lambdaFasta);

    // Open BAM file
    BamFile bamFile(singleInsertionBam);
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
    EXPECT_EQ("----------------------------------------------------AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
        r.ReferenceSubsequence(record, Orientation::NATIVE, true));
    EXPECT_EQ("AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
        r.ReferenceSubsequence(record, Orientation::NATIVE, true, true));
    EXPECT_EQ("----------------------------------------------------AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
        r.ReferenceSubsequence(record, Orientation::GENOMIC, true));
    EXPECT_EQ("AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
        r.ReferenceSubsequence(record, Orientation::GENOMIC, true, true));
    record = *it++;
    EXPECT_EQ("AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA----------------------------------------------------",
        r.ReferenceSubsequence(record, Orientation::GENOMIC, true));
    EXPECT_EQ("----------------------------------------------------TTGCCGCTGTT-ACCGTGCTGCGATCTTCTGCCATCGACGGACGTCCCACATTGGTGACTT",
        r.ReferenceSubsequence(record, Orientation::NATIVE, true));
    EXPECT_EQ("AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGT-AACAGCGGCAA",
        r.ReferenceSubsequence(record, Orientation::GENOMIC, true, true));
    EXPECT_EQ("TTGCCGCTGTT-ACCGTGCTGCGATCTTCTGCCATCGACGGACGTCCCACATTGGTGACTT",
        r.ReferenceSubsequence(record, Orientation::NATIVE, true, true));

    // {
    //     std::stringstream output;
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
    IndexedFastaReader r(lambdaFasta);

    EXPECT_TRUE(r.HasSequence("lambda_NEB3011"));
    EXPECT_FALSE(r.HasSequence("dog"));
    EXPECT_EQ(1, r.NumSequences());
    EXPECT_EQ(48502, r.SequenceLength("lambda_NEB3011"));

    string seq = r.Subsequence("lambda_NEB3011:0-10");
    EXPECT_EQ("GGGCGGCGAC", seq);

    string seq2 = r.Subsequence("lambda_NEB3011", 0, 10);
    EXPECT_EQ("GGGCGGCGAC", seq2);

    // subsequence extending beyond bounds returns clipped
    string seq3 = r.Subsequence("lambda_NEB3011", 48400, 48600);
    EXPECT_EQ(102, seq3.length());

    // bad subsequence
}

TEST(IndexedFastaReaderTests, Errors)
{
    IndexedFastaReader r(lambdaFasta);

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
    IndexedFastaReader r(lambdaFasta);
    std::vector<std::string> names = {"lambda_NEB3011"};

    // Test all-name request
    EXPECT_EQ(names, r.Names());

    // Test single-name query
    EXPECT_EQ(names[0], r.Name(0));

    // invalid name acces (out of range)
    EXPECT_THROW(r.Name(1), std::exception);
}
