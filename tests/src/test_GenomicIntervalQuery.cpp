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
#include <pbbam/GenomicIntervalQuery.h>
#include <iostream>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string inputBamFn = tests::Data_Dir + "/ex2.bam";

TEST(GenomicIntervalQueryTest, CountRecords)
{
    // open input BAM file
    BamFile bamFile(inputBamFn);
    EXPECT_TRUE(bamFile);

    const int id = bamFile.ReferenceId("seq1");
    EXPECT_TRUE(id != -1);

    // count records
    int count = 0;
    GenomicInterval interval(id, 0, 100);
    GenomicIntervalQuery query(interval, bamFile);
    EXPECT_TRUE(query);
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(39, count);

    // adjust interval and pass back in
    interval.Start(500);
    interval.Stop(600);
    query.Interval(interval);
    EXPECT_TRUE(query);
    count = 0;
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(166, count);

    // adjust again
    interval.Id(1);
    interval.Start(0);
    interval.Stop(100);
    query.Interval(interval);
    EXPECT_TRUE(query);
    count = 0;
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(83, count);

    // unknown ref (-1)
    interval.Id(-1);
    interval.Start(0);
    interval.Stop(100);
    query.Interval(interval);
    EXPECT_FALSE(query);                       // operator bool() returns false
    count = 0;
    for (const BamRecord& record : query) {    // iteration is still valid, just returns no data
        (void)record;
        ++count;
    }
    EXPECT_EQ(0, count);

    // adjust again - make sure we can read a real region after an invalid one
    interval.Id(1);
    interval.Start(0);
    interval.Stop(100);
    query.Interval(interval);
    EXPECT_TRUE(query);
    count = 0;
    for (const BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(83, count);

    // unknown ref (too high)
    interval.Id(2);
    interval.Start(0);
    interval.Stop(100);
    query.Interval(interval);
    EXPECT_FALSE(query);                       // operator bool() returns false
    count = 0;
    for (const BamRecord& record : query) {    // iteration is still valid, just returns no data
        (void)record;
        ++count;
    }
    EXPECT_EQ(0, count);
}

TEST(GenomicIntervalQueryTest, NonConstBamRecord)
{
    // open input BAM file
    BamFile bamFile(inputBamFn);
    EXPECT_TRUE(bamFile);

    const int id = bamFile.ReferenceId("seq1");
    EXPECT_TRUE(id != -1);

    // count records
    int count = 0;
    GenomicInterval interval(id, 0, 100);
    GenomicIntervalQuery query(interval, bamFile);
    EXPECT_TRUE(query);
    for (BamRecord& record : query) {
        (void)record;
        ++count;
    }
    EXPECT_EQ(39, count);
}

//TEST(GenomicIntervalQueryTest, WorksWithBamRecordImpl)
//{
//    // open input BAM file
//    BamFile bamFile(inputBamFn);
//    EXPECT_TRUE(bamFile);

//    const int id = bamFile.ReferenceId("seq1");
//    EXPECT_TRUE(id != -1);

//    // count records
//    int count = 0;
//    GenomicInterval interval(id, 0, 100);
//    GenomicIntervalQuery query(interval, bamFile);
//    EXPECT_TRUE(query);
//    for (const BamRecordImpl& record : query) {
//        (void)record;
//        ++count;
//    }
//    EXPECT_EQ(39, count);
//}

//TEST(GenomicIntervalQueryTest, WorksWithNonConstBamRecordImpl)
//{
//    // open input BAM file
//    BamFile bamFile(inputBamFn);
//    EXPECT_TRUE(bamFile);

//    const int id = bamFile.ReferenceId("seq1");
//    EXPECT_TRUE(id != -1);

//    // count records
//    int count = 0;
//    GenomicInterval interval(id, 0, 100);
//    GenomicIntervalQuery query(interval, bamFile);
//    EXPECT_TRUE(query);
//    for (BamRecordImpl& record : query) {
//        (void)record;
//        ++count;
//    }
//    EXPECT_EQ(39, count);
//}

// add special cases as needed
