// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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

#include <gtest/gtest.h>

#ifdef PBBAM_TESTING
#define private public
#endif

#include "TestData.h"

#include <pbbam/FastqReader.h>
#include <pbbam/FastqSequence.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

static void CheckSequence(const size_t index, const FastqSequence& seq)
{
    SCOPED_TRACE("checking Fastq seq:" + std::to_string(index));
    switch (index) {
        case 0 :
            EXPECT_EQ("1", seq.Name());
            EXPECT_EQ("TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
                      "ACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGA"
                      "GGAGAACGCAACTCCGCCGGCGCAGGCG", seq.Bases());
            EXPECT_EQ("[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                      "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                      "[[[[[[[[[[[[[[[[[[[[[[[[[[[[", seq.Qualities().Fastq());
            break;

        case 1 :
            EXPECT_EQ("2", seq.Name());
            EXPECT_EQ("TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
                      "ACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGA"
                      "GGAGAACGCAAC", seq.Bases());
            EXPECT_EQ("[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                      "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
                      "[[[[[[[[[[[[", seq.Qualities().Fastq());
            break;

        case 2 :
            EXPECT_EQ("3", seq.Name());
            EXPECT_EQ("TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA"
                      "ACCCTAACCCTAACACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA"
                      "ACCCTAACCCCTAACCCTAACCCT", seq.Bases());
            EXPECT_EQ("]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                      "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                      "]]]]]]]]]]]]]]]]]]]]]]]]", seq.Qualities().Fastq());
            break;

        default:
            ASSERT_TRUE(false); // invalid index
    }
}

TEST(FastqSequenceTest, BasicConstructorsOk)
{
    FastqSequence seq1{ "1", "GATTACA", "[[[[[[["};
    EXPECT_EQ("1",       seq1.Name());
    EXPECT_EQ("GATTACA", seq1.Bases());
    EXPECT_EQ("[[[[[[[", seq1.Qualities().Fastq());

    const auto quals = vector<uint8_t>{ 58,58,58,58,58,58,58 };
    FastqSequence seq2{ "1", "GATTACA", QualityValues{quals} };
    EXPECT_EQ("1",       seq2.Name());
    EXPECT_EQ("GATTACA", seq2.Bases());
    EXPECT_EQ("[[[[[[[", seq2.Qualities().Fastq());
}

TEST(FastqReaderTest, IterableOk)
{
    const string fn = tests::GeneratedData_Dir + "/normal.fq";
    FastqReader reader{ fn };

    size_t count = 0;
    FastqSequence seq;
    while (reader.GetNext(seq)) {
        CheckSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(FastqReaderTest, ReadAllOk)
{
    const string fn = tests::GeneratedData_Dir + "/normal.fq";

    size_t count = 0;
    for (const auto& seq : FastqReader::ReadAll(fn)) {
        CheckSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(3, count);
}

