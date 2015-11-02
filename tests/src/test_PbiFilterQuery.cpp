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
#include <pbbam/PbiFilterQuery.h>
#include <algorithm>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(PbiFilterQueryTest, QueryOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/test_group_query/test2.bam" } };

    {
        int count = 0;
        PbiFilterQuery query( PbiQueryLengthFilter{ 500, Compare::GREATER_THAN_EQUAL}, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_GE((r.QueryEnd() - r.QueryStart()), 500);
        }
        EXPECT_EQ(3, count);
    }
    {
        // all records aligned to reverse strand && pos >= 9200
        const auto filter = PbiFilter::Intersection(
        {
            PbiAlignedStrandFilter{Strand::REVERSE},
            PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}
        });

        int count = 0;
        PbiFilterQuery query(filter, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_EQ(Strand::REVERSE, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"), r.FullName());
        }
        EXPECT_EQ(1, count);
    }
    {
        // all records aligned to forward strand && pos >= 9200
        const auto filter = PbiFilter::Intersection(
        {
            PbiAlignedStrandFilter{Strand::FORWARD},
            PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}
        });

        int count = 0;
        PbiFilterQuery query(filter, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_EQ(Strand::FORWARD, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2114_2531"), r.FullName());
        }
        EXPECT_EQ(1, count);
    }
    {
        // all records from RG ("b89a4406") with numMatches >= 1200
        const auto filter = PbiFilter::Intersection(
        {
            PbiReadGroupFilter{"b89a4406"},
            PbiNumMatchesFilter{1200, Compare::GREATER_THAN_EQUAL}
        });

        int count = 0;
        PbiFilterQuery query(filter, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_EQ(string("b89a4406"), r.ReadGroupId());
            EXPECT_GE((r.NumMatches()), 1200);
            if (count == 1)
                EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055"), r.FullName());
            else if (count == 2)
                EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/4101_5571"), r.FullName());
        }
        EXPECT_EQ(2, count);
    }
}

TEST(PbiFilterQueryTest, ZmwRangeFromDatasetOk)
{
    const auto expectedMovieName = string{ "m150404_101626_42267_c100807920800000001823174110291514_s1_p0" };

    DataSet ds(tests::Data_Dir + "/chunking/chunking.subreadset.xml");
    EXPECT_EQ(3, ds.BamFiles().size());

    { // movie name

        int count = 0;
        PbiFilterQuery query{ PbiMovieNameFilter{expectedMovieName}, ds };
        for (const BamRecord& r : query) {
            EXPECT_EQ(expectedMovieName, r.MovieName());
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    
    { // sequencing chemistries
        set<string> chems{ ds.SequencingChemistries() };
        set<string> expected{ "P6-C4" };
        EXPECT_TRUE(equal(chems.begin(), chems.end(), expected.begin()));
    }

    { // min ZMW

        int count = 0;
        PbiFilterQuery query{ PbiZmwFilter{54, Compare::GREATER_THAN}, ds };
        for (const BamRecord& r : query) {
            EXPECT_GT(r.HoleNumber(), 54);
            ++count;
        }
        EXPECT_EQ(1220, count);
    }

    { // max ZMW

        int count = 0;
        PbiFilterQuery query{ PbiZmwFilter{1816, Compare::LESS_THAN}, ds };
        for (const BamRecord& r : query) {
            EXPECT_LT(r.HoleNumber(),1816);
            ++count;
        }
        EXPECT_EQ(150, count);
    }

    { // put all together, from DataSet XML

        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
        int count = 0;
        for (const BamRecord& r : query) {
            EXPECT_EQ(expectedMovieName, r.MovieName());
            const auto zmw = r.HoleNumber();
            EXPECT_GT(zmw, 54);
            EXPECT_LT(zmw, 1816);
            ++count;
        }
        EXPECT_EQ(150, count);
    }
    { // empty filter object - should return all records from the same dataset

        PbiFilterQuery query(PbiFilter{ }, ds);
        int count = 0;
        for (const BamRecord& r : query) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    { // no <Filters> element present at all

        const DataSet ds(tests::Data_Dir + "/chunking/chunking_missingfilters.subreadset.xml");
        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
        int count = 0;
        for (const BamRecord& r : query) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    { // <Filters> element contains no child <Filter> elements

        const DataSet ds(tests::Data_Dir + "/chunking/chunking_emptyfilters.subreadset.xml");
        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
        int count = 0;
        for (const BamRecord& r : query) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
}
