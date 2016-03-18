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
#include <pbbam/QNameQuery.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

static const string dataDir = tests::Data_Dir + "/group/";
static const string test1fn = string(dataDir) + "test1.bam";
static const string test2fn = string(dataDir) + "test2.bam";
static const string test3fn = string(dataDir) + "test3.bam";

static
void TestQNameQuery(const string& fn, const vector<int>& expected)
{
    EXPECT_NO_THROW(
    {
        vector<int> counts;
        QNameQuery qQuery(fn);
        for (const vector<BamRecord>& records : qQuery)
            counts.push_back(records.size());
        EXPECT_EQ(expected, counts);
    });
}

static
void TestNoneConstQNameQuery(const string& fn, const vector<int>& expected)
{
    EXPECT_NO_THROW(
    {
        vector<int> counts;
        QNameQuery qQuery(fn);
        for (vector<BamRecord>& records : qQuery)
            counts.push_back(records.size());
        EXPECT_EQ(expected, counts);
    });
}

TEST(QNameQueryTest, CountQSizes)
{
    // test case 1 has exactly one bamRecord.
    string fn = test1fn;
    vector<int> expected({1});
    TestQNameQuery(fn, expected);
    TestNoneConstQNameQuery(fn, expected);

    // test case 2 has bamRecords of four subreads.
    fn = test2fn;
    expected = {1, 1, 1, 1};
    TestQNameQuery(fn, expected);
    TestNoneConstQNameQuery(fn, expected);

    fn = test3fn;
    expected = {2,1,1,1,1,1,1,2,1,1,1};
    TestQNameQuery(fn, expected);
    TestNoneConstQNameQuery(fn, expected);
}

