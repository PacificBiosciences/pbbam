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
#include <pbbam/EntireFileQuery.h>
#include <pbbam/BamWriter.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string inputBamFn = tests::Data_Dir + "/ex2.bam";

TEST(EntireFileQueryTest, CountRecords)
{
    EXPECT_NO_THROW(
    {
        // open input BAM file
        BamFile bamFile(inputBamFn);

        // count records
        int count = 0;
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile) {
            (void)record;
            ++count;
        }

        EXPECT_EQ(3307, count);
    });
}

TEST(EntireFileQueryTest, NonConstBamRecord)
{
    EXPECT_NO_THROW(
    {
        // open input BAM file
        BamFile bamFile(inputBamFn);

        // count records
        int count = 0;
        EntireFileQuery entireFile(bamFile);
        for (BamRecord& record : entireFile) {
            (void)record;
            ++count;
        }

        EXPECT_EQ(3307, count);
    });
}

TEST(BamRecordTest, HandlesDeletionOK)
{
    // this file raised no error in Debug mode, but segfaulted when
    // trying to access the aligned qualities in Release mode

    EXPECT_NO_THROW(
    {
        // open input BAM file
        const string problemBamFn = tests::Data_Dir + "/segfault.bam";
        BamFile bamFile(problemBamFn);

        // count records
        int count = 0;
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile) {

            const auto rawQualities     = record.Qualities(Orientation::GENOMIC, false);
            const auto alignedQualities = record.Qualities(Orientation::GENOMIC, true);

            const string rawExpected =
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

            // 1=1D98=
            const string alignedExpected =
                "I!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

            EXPECT_EQ(rawExpected,     rawQualities.Fastq());
            EXPECT_EQ(alignedExpected, alignedQualities.Fastq());

            ++count;
        }

        EXPECT_EQ(1, count);
    });
}


TEST(BamRecordTest, ReferenceName)
{
    // check reference name of first record
//    {
        const string exampleBam  = tests::Data_Dir + "/ex2.bam";
        BamFile bamFile(exampleBam);
        EntireFileQuery records(bamFile);

        auto it = records.begin();
        auto record = *it;

//        EXPECT_EQ("seq1", records.begin()->ReferenceName());
//    }

//    // unmapped records have no reference name, should throw
//    {
//        const string exampleBam  = tests::Data_Dir + "/unmap1.bam";
//        BamFile bamFile(exampleBam);
//        staging::EntireFileQuery records(bamFile);

//        EXPECT_THROW(records.begin()->ReferenceName(), std::exception);
//    }
}

TEST(ZZZ_LastTest_Perhaps, CreateSmallerTestBam) 
{
    BamFile file("/mnt/secondary-siv/mdsmith/simulationStudies/lambdaToMerge/10Gb/pbalchemysim.all.pbalign.bam");
    BamWriter writer("smaller.bam", file.Header());
    EntireFileQuery query(file);
    int count = 0;
    for (const BamRecord& b : query) {
        writer.Write(b);
        ++count;
        if (count == 25000)
            break;
    } 
}

// add add'l special cases as needed
