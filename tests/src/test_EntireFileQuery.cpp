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

#include <pbbam/EntireFileQuery.h>
#include <pbbam/BamWriter.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace EntireFileQueryTests {

const string inputBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";

} // namespace EntireFileQueryTests

TEST(EntireFileQueryTest, CountRecords)
{
    EXPECT_NO_THROW(
    {
        BamFile bamFile(EntireFileQueryTests::inputBamFn);
        int count = 0;
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile) {
            UNUSED(record);
            ++count;
        }

        EXPECT_EQ(4, count);
    });
}

TEST(EntireFileQueryTest, NonConstBamRecord)
{
    EXPECT_NO_THROW(
    {
        BamFile bamFile(EntireFileQueryTests::inputBamFn);
        int count = 0;
        EntireFileQuery entireFile(bamFile);
        for (BamRecord& record : entireFile) {
            UNUSED(record);
            ++count;
        }

        EXPECT_EQ(4, count);
    });
}

TEST(BamRecordTest, HandlesDeletionOK)
{
    // this file raised no error in Debug mode, but segfaulted when
    // trying to access the aligned qualities in Release mode

    const string problemBamFn = PbbamTestsConfig::Data_Dir + "/segfault.bam";
    BamFile bamFile(problemBamFn);
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
}


TEST(BamRecordTest, ReferenceName)
{
    {   // check reference name of first record
        const string exampleBam  = PbbamTestsConfig::Data_Dir + "/aligned.bam";
        BamFile bamFile(exampleBam);
        EntireFileQuery records(bamFile);
        auto firstIter = records.begin();
        auto& firstRecord = *firstIter;
        ASSERT_TRUE(firstRecord.IsMapped());
        EXPECT_EQ("lambda_NEB3011", firstRecord.ReferenceName());
    }

    {   // unmapped records have no reference name, should throw
        const string exampleBam  = PbbamTestsConfig::Data_Dir + "/unmap1.bam";
        BamFile bamFile(exampleBam);
        EntireFileQuery records(bamFile);
        auto firstIter = records.begin();
        auto& firstRecord = *firstIter;
        ASSERT_FALSE(firstRecord.IsMapped());
        EXPECT_THROW(firstRecord.ReferenceName(), std::runtime_error);
    }
}
