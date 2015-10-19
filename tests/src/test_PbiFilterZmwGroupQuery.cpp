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
#include <pbbam/PbiFilterZmwGroupQuery.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

static const string testChunking = tests::Data_Dir + "/chunking/chunking.subreadset.xml";
static const string testNoFilter = tests::Data_Dir + "/chunking/nofilter.subreadset.xml";

static
void TestPbiFilterZmwGroupQuery(const string& fn, 
                                const vector<size_t>& expected,
                                const uint32_t min_zmw,
                                const uint32_t max_zmw)
{
    EXPECT_NO_THROW(
    {
        PbiFilterZmwGroupQuery qQuery(fn);

        vector<size_t> counts;
        for (const vector<BamRecord>& records : qQuery)
        {
            counts.push_back(records.size());
            EXPECT_GT(records.size(), 0);
            string movieName = records[0].MovieName();
            uint32_t holeNumber = records[0].HoleNumber();
            EXPECT_TRUE(holeNumber >= min_zmw);
            EXPECT_TRUE(holeNumber <= max_zmw);
            for (const BamRecord & record: records) {
                EXPECT_EQ(holeNumber, record.HoleNumber());
                EXPECT_EQ(movieName, record.MovieName());
            }
        }
        EXPECT_EQ(expected, counts);
    });

}

static
void TestNoneConstPbiFilterZmwGroupQuery(const string& fn,
                                         const vector<size_t>& expected,
                                         const uint32_t min_zmw,
                                         const uint32_t max_zmw)
{
    EXPECT_NO_THROW(
    {
        PbiFilterZmwGroupQuery qQuery(fn);

        vector<size_t> counts;
        for (vector<BamRecord>& records : qQuery) {
            counts.push_back(records.size());
            EXPECT_GT(records.size(), 0);
            string movieName = records[0].MovieName();
            uint32_t holeNumber = records[0].HoleNumber();
            EXPECT_TRUE(holeNumber >= min_zmw);
            EXPECT_TRUE(holeNumber <= max_zmw);
            for (BamRecord & record: records) {
                EXPECT_EQ(holeNumber, record.HoleNumber());
                EXPECT_EQ(movieName, record.MovieName());
            }
        }
        EXPECT_EQ(expected, counts);
    });

}

TEST(PbiFilterZmwGroupQueryTest, GetNext)
{
    string fn = testChunking;
    vector<size_t> expected({2, 21, 13, 1, 5, 13, 1, 34, 12, 2, 20, 5, 3, 7, 11});
    const uint32_t min_zmw = 55;
    const uint32_t max_zmw = 1816;
    TestPbiFilterZmwGroupQuery(fn, expected, min_zmw, max_zmw);
    TestNoneConstPbiFilterZmwGroupQuery(fn, expected, min_zmw, max_zmw);
}

TEST(PbiFilterZmwGroupQueryTest, NoFilter)
{
    string fn = testNoFilter;
    vector<size_t> expected({2,21,13,1,5,13,1,34,12,2,20,5,3,7,11,14,6,8,23,53,17,21,7,5,35,3,26,6,21,37,26,59,2,6,30,34,32,2,14,3,24,1,15,1,12,26,6,3,1,9,3,21,12,10,24,3,6,1,6,17,34,11,24,4,11,1,10,8,10,20,3,4,6,27,5,2,21,3,14,1,9,5,30,37,6,1,26,7,7,32});
    const uint32_t min_zmw = 0;
    const uint32_t max_zmw = 1000000;
    TestPbiFilterZmwGroupQuery(fn, expected, min_zmw, max_zmw);
    TestNoneConstPbiFilterZmwGroupQuery(fn, expected, min_zmw, max_zmw);
}
