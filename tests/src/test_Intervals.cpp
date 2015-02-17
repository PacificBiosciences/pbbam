// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#include <gtest/gtest.h>
#include <pbbam/GenomicInterval.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(IntervalTest, Constructors)
{
    Interval<uint32_t> empty;
    Interval<uint32_t> singleton(4);
    Interval<uint32_t> normal(5, 8);

    EXPECT_EQ(0, empty.Start());
    EXPECT_EQ(0, empty.Stop());

    EXPECT_EQ(4, singleton.Start());
    EXPECT_EQ(5, singleton.Stop());

    EXPECT_EQ(5, normal.Start());
    EXPECT_EQ(8, normal.Stop());

    // TODO: check out-of-order intervals, etc
}

TEST(IntervalTest, EqualityTest)
{
    Interval<uint32_t> empty;
    Interval<uint32_t> empty2;

    Interval<uint32_t> singleton(4);
    Interval<uint32_t> sameAsSingleton(4,5);

    Interval<uint32_t> normal(5, 8);
    Interval<uint32_t> sameAsNormal(5, 8);

    Interval<uint32_t> different(20, 40);

    // self-equality
    EXPECT_TRUE(empty == empty);
    EXPECT_TRUE(singleton == singleton);
    EXPECT_TRUE(normal == normal);
    EXPECT_TRUE(different == different);

    // same values equality
    EXPECT_TRUE(empty == empty2);
    EXPECT_TRUE(singleton == sameAsSingleton);
    EXPECT_TRUE(normal == sameAsNormal);

    // different values
    EXPECT_FALSE(empty == singleton);
    EXPECT_FALSE(empty == normal);
    EXPECT_FALSE(empty == different);
    EXPECT_FALSE(singleton == normal);
    EXPECT_FALSE(normal == different);
}

TEST(IntervalTest, Copy)
{
    Interval<uint32_t> interval1(5, 8);
    Interval<uint32_t> interval2(interval1);
    Interval<uint32_t> interval3 = interval1;

    EXPECT_TRUE(interval1 == interval1);
    EXPECT_TRUE(interval1 == interval2);
    EXPECT_TRUE(interval1 == interval3);
}

TEST(IntervalTest, Modifier)
{
    Interval<uint32_t> interval1(5, 8);
    Interval<uint32_t> interval2(interval1);
    interval2.Start(2);
    interval2.Stop(10);

    EXPECT_FALSE(interval1 == interval2);
    EXPECT_EQ(2, interval2.Start());
    EXPECT_EQ(10, interval2.Stop());
}

TEST(IntervalTest, CoverTest)
{
    Interval<uint32_t> interval1(2, 4);
    Interval<uint32_t> interval2(3, 5);
    Interval<uint32_t> interval3(6, 8);
    Interval<uint32_t> interval4(1, 7);
    Interval<uint32_t> interval5(5, 8);

    EXPECT_TRUE(interval1.Covers(interval1));    // self-cover: a.covers(a)
    EXPECT_TRUE(interval1.CoveredBy(interval1)); // self-cover: a.coveredBy(a)

    EXPECT_TRUE(interval2.CoveredBy(interval4)); // a.coveredBy(b)
    EXPECT_TRUE(interval4.Covers(interval2));    // thus b.covers(a)
    EXPECT_FALSE(interval2 == interval4);        // if a != b
    EXPECT_FALSE(interval2.Covers(interval4));   // then !a.covers(b)

    EXPECT_FALSE(interval2.Covers(interval3));    // completely disjoint
    EXPECT_FALSE(interval3.Covers(interval2));
    EXPECT_FALSE(interval2.CoveredBy(interval3));
    EXPECT_FALSE(interval3.CoveredBy(interval2));

    EXPECT_FALSE(interval2.Covers(interval5));    // a.stop == b.start
    EXPECT_FALSE(interval2.CoveredBy(interval5));

    EXPECT_TRUE(interval5.Covers(interval3));    // shared endpoint, start contained, thus a.covers(b)
    EXPECT_TRUE(interval3.CoveredBy(interval5)); // and b.coveredBy(a)
}

TEST(IntervalTest, IntersectTest)
{
    Interval<uint32_t> interval1(2, 4);
    Interval<uint32_t> interval2(3, 5);
    Interval<uint32_t> interval3(6, 8);
    Interval<uint32_t> interval4(1, 7);
    Interval<uint32_t> interval5(5, 8);

    EXPECT_TRUE(interval1.Intersects(interval1)); // self-intersection: a.intersects(a)

    EXPECT_TRUE(interval1.Intersects(interval2)); // if a.intersects(b)
    EXPECT_TRUE(interval2.Intersects(interval1)); // then b.intersects(a)

    EXPECT_TRUE(interval4.Covers(interval1));      // if b.covers(a),
    EXPECT_TRUE(interval1.Intersects(interval4));  // then a.intersects(b)
    EXPECT_TRUE(interval4.Intersects(interval1));  // and b.intersects(a)

    EXPECT_FALSE(interval2.Intersects(interval3)); // b.start > a.stop (obvious disjoint)
    EXPECT_FALSE(interval2.Intersects(interval5)); // b.start == a.stop (intervals are right open, so disjoint)
}

TEST(IntervalTest, LengthTest)
{
    Interval<uint32_t> interval1(2, 4);
    Interval<uint32_t> interval2(3, 5);
    Interval<uint32_t> interval3(6, 8);
    Interval<uint32_t> interval4(1, 7);
    Interval<uint32_t> interval5(5, 8);

    EXPECT_EQ(2, interval1.Length());
    EXPECT_EQ(2, interval2.Length());
    EXPECT_EQ(2, interval3.Length());
    EXPECT_EQ(6, interval4.Length());
    EXPECT_EQ(3, interval5.Length());

    // TODO: check out-of-order intervals, etc
}

TEST(GenomicIntervalTest, Constructors)
{
    GenomicInterval empty;
    GenomicInterval refAndCoords(0, 100, 200);

    EXPECT_EQ(-1, empty.Id());
    EXPECT_EQ(0,  empty.Start());
    EXPECT_EQ(0,  empty.Stop());

    EXPECT_EQ(0,   refAndCoords.Id());
    EXPECT_EQ(100, refAndCoords.Start());
    EXPECT_EQ(200, refAndCoords.Stop());
}

TEST(GenomicIntervalTest, Copy)
{
    GenomicInterval interval1(1, 10, 20);
    GenomicInterval interval2(interval1);
    GenomicInterval interval3 = interval1;

    EXPECT_TRUE(interval1 == interval1);
    EXPECT_TRUE(interval1 == interval2);
    EXPECT_TRUE(interval1 == interval3);
}

TEST(GenomicIntervalTest, Modifiers)
{
    GenomicInterval interval1(1, 10, 20);
    GenomicInterval interval2(interval1);
    interval2.Id(5);
    interval2.Start(2);
    interval2.Stop(10);

    EXPECT_FALSE(interval1 == interval2);
    EXPECT_EQ(5,  interval2.Id());
    EXPECT_EQ(2,  interval2.Start());
    EXPECT_EQ(10, interval2.Stop());
}

TEST(GenomicIntervalTest, CoverTest)
{
    GenomicInterval interval1(0, 2, 4);
    GenomicInterval interval2(0, 3, 5);
    GenomicInterval interval3(0, 6, 8);
    GenomicInterval interval4(0, 1, 7);
    GenomicInterval interval5(0, 5, 8);

    // same as interval2, but different ref
    GenomicInterval interval6(interval2);
    interval6.Id(1);

    EXPECT_TRUE(interval1.Covers(interval1));    // self-cover: a.covers(a)
    EXPECT_TRUE(interval1.CoveredBy(interval1)); // self-cover: a.coveredBy(a)

    EXPECT_TRUE(interval2.CoveredBy(interval4)); // a.coveredBy(b)
    EXPECT_TRUE(interval4.Covers(interval2));    // thus b.covers(a)
    EXPECT_FALSE(interval2 == interval4);        // if a != b
    EXPECT_FALSE(interval2.Covers(interval4));   // then !a.covers(b)

    EXPECT_FALSE(interval6.CoveredBy(interval4)); // interval 6 has same start/stop as 2, w/ different ref
    EXPECT_FALSE(interval4.Covers(interval6));    //
    EXPECT_FALSE(interval6 == interval4);         //
    EXPECT_FALSE(interval6.Covers(interval4));    //

    EXPECT_FALSE(interval2.Covers(interval3));    // completely disjoint
    EXPECT_FALSE(interval3.Covers(interval2));
    EXPECT_FALSE(interval2.CoveredBy(interval3));
    EXPECT_FALSE(interval3.CoveredBy(interval2));

    EXPECT_FALSE(interval2.Covers(interval5));    // a.stop == b.start
    EXPECT_FALSE(interval2.CoveredBy(interval5));

    EXPECT_TRUE(interval5.Covers(interval3));    // shared endpoint, start contained, thus a.covers(b)
    EXPECT_TRUE(interval3.CoveredBy(interval5)); // and b.coveredBy(a)
}

