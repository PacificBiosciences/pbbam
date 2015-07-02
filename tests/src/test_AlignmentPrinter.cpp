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

#include "pbbam/AlignmentPrinter.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/IndexedFastaReader.h"

#include <iostream>
#include <sstream>
#include <string>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string lambdaFasta = tests::Data_Dir + "/lambdaNEB.fa";
const string singleInsertionBam = tests::Data_Dir + "/aligned.bam";

TEST(AlignmentPrinterTest, Print)
{
    IndexedFastaReader r(lambdaFasta);

    BamFile bamFile(singleInsertionBam);
    EntireFileQuery bamQuery(bamFile);

    auto it = bamQuery.begin();
    

    // std::cerr << record.AlignedStart() << std::endl;
    // std::cerr << record.Sequence(Orientation::GENOMIC, true) << std::endl;
    // std::cerr << record.Sequence(Orientation::GENOMIC, true, true) << std::endl;

    AlignmentPrinter pretty(r);

    // std::string expected = 
    // "Read        : singleInsertion2\n"
    // "Reference   : lambda_NEB3011\n"
    // "\n"
    // "Read-length : 49\n"
    // "Concordance : 0.96\n"
    // "\n"
    // "   GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGGACTGGCTGAT\n"
    // "   |||||||| ||||||||||||||||||| |||||||||||||||||||||\n"
    // "   GGCTGCAG-GTACAGCGGTCAGGAGGCCAATTGATGCCGGACTGGCTGAT\n";
    // EXPECT_EQ(expected, pretty.Print(record, Orientation::NATIVE));

    auto record = *it++;
    std::cerr << pretty.Print(record, Orientation::GENOMIC);
    std::cerr << std::endl << std::endl;
    record = *it++;
    std::cerr << pretty.Print(record, Orientation::GENOMIC);
    std::cerr << std::endl << std::endl;
    record = *it++;
    std::cerr << pretty.Print(record, Orientation::GENOMIC);
    std::cerr << std::endl << std::endl;
    record = *it++;
    std::cerr << pretty.Print(record, Orientation::GENOMIC);
    std::cerr << std::endl << std::endl;
}
