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
#include <string>

#include <gtest/gtest.h>

#define private public

#include "PbbamTestData.h"

#include <pbbam/GenomicIntervalQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace GenomicIntervalQueryTests {
    const string inputBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
} // namespace GenomicIntervalQueryTests

TEST(GenomicIntervalQueryTest, ReuseQueryAndCountRecords)
{
    const string rname = "lambda_NEB3011";

    BamFile bamFile(GenomicIntervalQueryTests::inputBamFn);

    // setup with normal interval
    int count = 0;
    GenomicInterval interval(rname, 5000, 6000);
    GenomicIntervalQuery query(interval, bamFile);
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(2, count);

    // adjust interval and pass back in
    count = 0;
    interval.Start(9300);
    interval.Stop(9400);
    query.Interval(interval);
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(2, count);

    // adjust again (empty region)
    count = 0;
    interval.Name(rname);
    interval.Start(1000);
    interval.Stop(2000);
    query.Interval(interval);
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(0, count);

    // unknown ref
    count = 0;
    interval.Name("does not exist");
    interval.Start(0);
    interval.Stop(100);
    EXPECT_THROW(query.Interval(interval), std::runtime_error);
    for (const BamRecord& record : query) { // iteration is still safe, just returns no data
        (void)record;
        ++count;
    }
    EXPECT_EQ(0, count);

    // adjust again - make sure we can read a real region after an invalid one
    interval.Name(rname);
    interval.Start(5000);
    interval.Stop(6000);
    query.Interval(interval);
    count = 0;
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(2, count);
}

TEST(GenomicIntervalQueryTest, NonConstBamRecord)
{
    EXPECT_NO_THROW(
    {
        BamFile bamFile(GenomicIntervalQueryTests::inputBamFn);
        int count = 0;

        GenomicInterval interval("lambda_NEB3011", 8000, 10000);
        GenomicIntervalQuery query(interval, bamFile);
        for (BamRecord& record : query) {
            (void)record;
            ++count;
        }
        EXPECT_EQ(2, count);
    });
}

TEST(GenomicIntervalQueryTest,  MissingBaiShouldThrow)
{
    GenomicInterval interval("lambda_NEB3011", 0, 100);
    const string phi29Bam = PbbamTestsConfig::Data_Dir + "/phi29.bam";
    const string hasBaiBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";

    {   // single file, missing BAI
        EXPECT_THROW(GenomicIntervalQuery query(interval, phi29Bam), std::runtime_error);
    }

    {   // from dataset, all missing BAI
        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        EXPECT_THROW(GenomicIntervalQuery query(interval, ds), std::runtime_error);
    }

    {   // from dataset, mixed BAI presence
        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(ExternalResource("PacBio.AlignmentFile.AlignmentBamFile", hasBaiBam));
        EXPECT_THROW(GenomicIntervalQuery query(interval, ds), std::runtime_error);
    }
}
