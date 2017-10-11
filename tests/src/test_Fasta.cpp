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

#include <cstddef>

#include <gtest/gtest.h>

#define private public

#include "PbbamTestData.h"

#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace FastaTests {

static void CheckSequence(const size_t index, const FastaSequence& seq)
{
    SCOPED_TRACE("checking FASTA seq:" + std::to_string(index));
    switch (index) {
        case 0 :
            EXPECT_EQ("1", seq.Name());
            EXPECT_EQ("TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGGCGCAGGCG", seq.Bases());
            break;

        case 1 :
            EXPECT_EQ("2", seq.Name());
            EXPECT_EQ("TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAAC", seq.Bases());
            break;

        case 2 :
            EXPECT_EQ("3", seq.Name());
            EXPECT_EQ("TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCT", seq.Bases());
            break;

        default:
            ASSERT_TRUE(false); // invalid index
    }
}

} // namespace FastaTests

TEST(FastaSequenceTest, BasicConstructorOk)
{
    FastaSequence seq{ "1", "GATTACA" };
    EXPECT_EQ("1",       seq.Name());
    EXPECT_EQ("GATTACA", seq.Bases());
}

TEST(FastaReaderTest, IterableOk)
{
    const string fn = PbbamTestsConfig::GeneratedData_Dir + "/normal.fa";
    FastaReader reader{ fn };

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
    const string fn = PbbamTestsConfig::GeneratedData_Dir + "/normal.fa";

    size_t count = 0;
    for (const auto& seq : FastaReader::ReadAll(fn)) {
        FastaTests::CheckSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(FastaSequenceQueryTest, FromFastaFilename)
{
    const string fn = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";

    {
        size_t count = 0;
        FastaSequenceQuery query{ fn } ;
        for (const auto& seq : query) {
            UNUSED(seq);
            ++count;
        }
        EXPECT_EQ(1, count);
    }

    {
        FastaSequenceQuery query{ fn };
        const auto first = query.cbegin();
        const auto& seq = *first;
        EXPECT_EQ("lambda_NEB3011", seq.Name());
    }

}

TEST(FastaSequenceQueryTest, FromDataSet)
{
    const string fn = PbbamTestsConfig::Data_Dir + "/referenceset.xml";

    {
        size_t count = 0;
        FastaSequenceQuery query{ fn } ;
        for (const auto& seq : query) {
            UNUSED(seq);
            ++count;
        }
        EXPECT_EQ(5, count);    // 1 from lambda, 4 from chimera
    }
    {
        FastaSequenceQuery query{ fn };
        const auto first = query.cbegin();
        const auto& seq = *first;
        EXPECT_EQ("lambda_NEB3011", seq.Name());
    }
}

