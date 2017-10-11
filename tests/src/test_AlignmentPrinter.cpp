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

#include <string>

#include <gtest/gtest.h>

#define private public

#include "PbbamTestData.h"

#include <pbbam/AlignmentPrinter.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/IndexedFastaReader.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace AlignmentPrinterTests {

const string lambdaFasta = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";
const string singleInsertionBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";

} // namespace AlignmentPrinterTests

TEST(AlignmentPrinterTest, Print)
{
    IndexedFastaReader r(AlignmentPrinterTests::lambdaFasta);
    AlignmentPrinter pretty(r);

    BamFile bamFile(AlignmentPrinterTests::singleInsertionBam);
    EntireFileQuery bamQuery(bamFile);
    auto it = bamQuery.begin();

    // funky formatting used to format alignments
    auto expected = string
    {
        "Read        : singleInsertion/100/0_49\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 49\n"
        "Concordance : 0.96\n"
        "\n"
        "5210 : GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGG : 5249\n"
        "       \x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| |\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| ||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "   0 : GGCTGCAG-GTACAGCGGTCAGGAGGCCAATTGATGCCGG :   39\n"
        "\n"
        "5249 : ACTGGCTGAT : 5259\n"
        "       |\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "  39 : ACTGGCTGAT :   49\n"
        "\n"
    };

    auto record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));

    expected = string 
    {
        "Read        : singleInsertion/200/0_49\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 49\n"
        "Concordance : 0.96\n"
        "\n"
        "5210 : GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGG : 5249\n"
        "       \x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| |\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| ||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "   0 : GGCTGCAG-GTACAGCGGTCAGGAGGCCAATTGATGCCGG :   39\n"
        "\n"
        "5249 : ACTGGCTGAT : 5259\n"
        "       |\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "  39 : ACTGGCTGAT :   49\n"
        "\n"
    };

    record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));

    expected = string 
    {
        "Read        : singleInsertion/100/0_111\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 59\n"
        "Concordance : 0.951\n"
        "\n"
        "9377 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCG : 9417\n"
        "       |||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||  |\n"
        "   0 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGA--G :   38\n"
        "\n"
        "9417 : CAGCACGGT-AACAGCGGCAA : 9437\n"
        "       |||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||| ||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||\n"
        "  38 : CAGCACGGTAAACAGCGGCAA :   59\n"
        "\n"
    };

    record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));

    expected = string 
    {
        "Read        : singleInsertion/100/0_111\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 59\n"
        "Concordance : 0.951\n"
        "\n"
        "9377 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCG : 9417\n"
        "       |||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||  |\n"
        "   0 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGA--G :   38\n"
        "\n"
        "9417 : CAGCACGGT-AACAGCGGCAA : 9437\n"
        "       |||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||| ||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||\n"
        "  38 : CAGCACGGTAAACAGCGGCAA :   59\n"
        "\n"
    };

    record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));
}
