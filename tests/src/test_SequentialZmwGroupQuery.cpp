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

// Author: Yuan Li

#include "TestData.h"
#include <gtest/gtest.h>
#include <pbbam/SequentialZmwGroupQuery.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

static const string dataDir = tests::Data_Dir + "/test_group_query/";
static const string testbam = string(dataDir) + "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.1.subreads.bam";
static const string testChunking = tests::Data_Dir + "/chunking/chunking.subreadset.xml";

static
void TestSequentialZmwGroupQuery(const string& fn, 
                                 const vector<int>& expected)
{
    EXPECT_NO_THROW(
    {
        vector<int> counts;
        SequentialZmwGroupQuery qQuery(fn);
        for (const vector<BamRecord>& records : qQuery)
        {
            counts.push_back(records.size());
            EXPECT_GT(records.size(), 0);
            string movieName = records[0].MovieName();
            uint32_t holeNumber = records[0].HoleNumber();
            for (const BamRecord & record: records) {
                EXPECT_EQ(holeNumber, record.HoleNumber());
                EXPECT_EQ(movieName, record.MovieName());
            }
        }
        EXPECT_EQ(expected, counts);
    });
}

static
void TestNoneConstSequentialZmwGroupQuery(const string& fn, 
                                          const vector<int>& expected)
{
    EXPECT_NO_THROW(
    {
        vector<int> counts;
        SequentialZmwGroupQuery qQuery(fn);
        for (vector<BamRecord>& records : qQuery) {
            counts.push_back(records.size());
            EXPECT_GT(records.size(), 0);
            string movieName = records[0].MovieName();
            uint32_t holeNumber = records[0].HoleNumber();
            for (BamRecord & record: records) {
                EXPECT_EQ(holeNumber, record.HoleNumber());
                EXPECT_EQ(movieName, record.MovieName());
            }
        }
        EXPECT_EQ(expected, counts);
    });
}

TEST(SequentialZmwGroupQueryTest, CountQSizes)
{
    string fn = testbam;
    vector<int> expected({2, 2, 10, 2, 3, 1, 2, 2, 3, 4, 1, 3, 1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 3, 8, 1, 3, 2, 1, 15, 2, 1, 3, 1, 2, 2, 1, 3, 3, 2, 2, 1, 2, 2, 1, 1, 1});
    TestSequentialZmwGroupQuery(fn, expected);
    TestNoneConstSequentialZmwGroupQuery(fn, expected);
} 

TEST(SequentialZmwGroupQueryTest, Chunking)
{
    string fn = testChunking;
    vector<int> expected({2,21,13,1,5,13,1,34,12,2,20,5,3,7,11,14,6,8,23,53,17,21,7,5,35,3,26,6,21,37,26,59,2,6,30,34,32,2,14,3,24,1,15,1,12,26,6,3,1,9,3,21,12,10,24,3,6,1,6,17,34,11,24,4,11,1,10,8,10,20,3,4,6,27,5,2,21,3,14,1,9,5,30,37,6,1,26,7,7,32});
    TestSequentialZmwGroupQuery(fn, expected);
}
