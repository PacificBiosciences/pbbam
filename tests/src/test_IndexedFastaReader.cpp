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

#ifdef PBBAM_TESTING
#define private public
#endif

#include "TestData.h"
#include <gtest/gtest.h>
#include <pbbam/IndexedFastaReader.h>
#include <iostream>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string lambdaFasta = tests::Data_Dir + "/lambdaNEB.fa";

TEST(IndexedFastaReaderTests, ReadLambda)
{
    IndexedFastaReader r;
    r.Open(lambdaFasta);

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
    IndexedFastaReader r;

    //
    // attempt access without "opening"
    //
    EXPECT_THROW(r.NumSequences(), std::exception);
    EXPECT_THROW(r.HasSequence("lambda_NEB3011"), std::exception);
    EXPECT_THROW(r.SequenceLength("lambda_NEB3011"), std::exception);
    EXPECT_THROW(r.Subsequence("lambda_NEB3011:0-10"), std::exception);

    //
    // invalid accesses after opening
    //
    r.Open(lambdaFasta);
    EXPECT_THROW(r.SequenceLength("dog"), std::exception);
    EXPECT_THROW(r.Subsequence("dog:0-10"), std::exception);
}







