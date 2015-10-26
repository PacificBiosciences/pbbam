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
#include <pbbam/PbiFilter.h>
#include <string>
#include <cstdio>
#include <cstdlib>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace tests {

// helper structs & methods

static
PbiRawData test2Bam_RawIndex(void)
{
    PbiRawData index;
    index.NumReads(4);

    PbiRawBasicData& subreadData = index.BasicData();
    subreadData.rgId_       = { -1197849594, -1197849594, -1197849594, -1197849594 };
    subreadData.qStart_     = { 2114, 2579, 4101, 5615 };
    subreadData.qEnd_       = { 2531, 4055, 5571, 6237 };
    subreadData.holeNumber_ = { 14743, 14743, 14743, 14743 };
    subreadData.readQual_   = { 0.901, 0.601, 0.901, 0.601 };
    subreadData.ctxtFlag_   = { 0, 0, 0, 0 };
    subreadData.fileOffset_ = { 35651584, 35655125, 35667128, 35679170 };

    PbiRawMappedData& mappedData = index.mappedData_;
    mappedData.tId_       = { 0, 0, 0, 0 };
    mappedData.tStart_    = { 9507, 8453, 8455, 9291 };
    mappedData.tEnd_      = { 9903, 9902, 9893, 9900 };
    mappedData.aStart_    = { 2130, 2581, 4102, 5619 };
    mappedData.aEnd_      = { 2531, 4055, 5560, 6237 };
    mappedData.revStrand_ = { 0, 1, 0, 1 };
    mappedData.mapQV_     = { 254, 254, 254, 254 };
    mappedData.nM_        = { 384, 1411, 1393, 598 };
    mappedData.nMM_       = { 0, 0, 0, 0 };

    PbiRawBarcodeData& barcodeData = index.barcodeData_;
    barcodeData.bcForward_ = { 0, 17, 256, 17 };
    barcodeData.bcReverse_ = { 1, 18, 257, 18 };
    barcodeData.bcQual_    = { 42, 80, 42, 110 };

    PbiRawReferenceData& referenceData = index.referenceData_;
    referenceData.entries_.emplace_back( 0, 0, 3 );
    referenceData.entries_.emplace_back( 1 );
    referenceData.entries_.emplace_back( PbiReferenceEntry::UNMAPPED_ID );

    return index;
}

static
PbiIndex createTestIndex(void)
{
    PbiIndex idx;
    idx.d_.reset(new ::PacBio::BAM::internal::PbiIndexPrivate(test2Bam_RawIndex()));
    return idx;
}

static const PbiIndex shared_index = createTestIndex();

static
void checkFilterInternals(const PbiFilter& filter,
                          const PbiFilter::CompositionType expectedType,
                          const size_t expectedNumChildren,
                          const size_t expectedLookupSize)
{
    EXPECT_EQ(expectedType,        filter.d_->type_);
    EXPECT_EQ(expectedNumChildren, filter.d_->filters_.size());
    EXPECT_EQ(expectedLookupSize,  filter.Lookup(tests::shared_index).size());
}

struct SimpleFilter
{
    IndexList Lookup(const PbiIndex&) const
    { return IndexList{ }; }
};

struct NoncompliantFilter { };

struct SortUniqueTestFilter
{
    IndexList Lookup(const PbiIndex&) const
    { return IndexList{ 2, 7, 0, 3, 4, 1, 8 }; }
};

struct SortUniqueTestFilter2
{
    IndexList Lookup(const PbiIndex&) const
    { return IndexList{ 3, 7, 5 }; }
};

static inline
PbiFilter emptyFilter(void)
{ return PbiFilter{ }; }

static inline
PbiFilter simpleFilter(void)
{ return PbiFilter{ SimpleFilter{ } }; }

} // namespace tests
} // namespace BAM
} // namespace PacBio

TEST(PbiFilterTest, DefaultCtorOk)
{
    auto filter = PbiFilter{ };
    tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 0, 4);
}

TEST(PbiFilterTest, CompositionOk)
{
    auto filter = PbiFilter{ };
    filter.Add(PbiFilter{ });
    tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, 4);
}

TEST(PbiFilterTest, CustomFilterOk)
{
    { // ctor
        auto filter = PbiFilter{ tests::SimpleFilter{ } };
        tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, 0);
    }
    { // Add
        auto filter = PbiFilter{ };
        filter.Add(tests::SimpleFilter{ });
        tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, 0);
    }

//    PbiFilter shouldNotCompile = PbiFilter{ tests::NoncompliantFilter{ } };                       // <-- when uncommented, should not compile
//    PbiFilter shouldNotCompileEither; shouldNotCompileEither.Add(tests::NoncompliantFilter{ });   // <-- when uncommented, should not compile
}

TEST(PbiFilterTest, CopyOk)
{
    { // empty
        const auto original = PbiFilter{ };

        PbiFilter copyCtor(original);
        PbiFilter copyAssign;
        copyAssign = original;

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 0, 4);
        tests::checkFilterInternals(copyCtor,   PbiFilter::INTERSECT, 0, 4);
        tests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 0, 4);

        EXPECT_EQ(4, original.Lookup(tests::shared_index).size());
        EXPECT_EQ(4, copyCtor.Lookup(tests::shared_index).size());
        EXPECT_EQ(4, copyAssign.Lookup(tests::shared_index).size());
    }
    { // with children
        const auto original = PbiFilter{ tests::SimpleFilter{ } };

        PbiFilter copyCtor(original);
        PbiFilter copyAssign;
        copyAssign = original;

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 1, 0);
        tests::checkFilterInternals(copyCtor,   PbiFilter::INTERSECT, 1, 0);
        tests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 1, 0);

        EXPECT_EQ(0, original.Lookup(tests::shared_index).size());
        EXPECT_EQ(0, copyCtor.Lookup(tests::shared_index).size());
        EXPECT_EQ(0, copyAssign.Lookup(tests::shared_index).size());
    }
}

TEST(PbiFilterTest, MoveOk)
{
    { // empty
        const auto original = tests::emptyFilter();

        PbiFilter moveCtor(tests::emptyFilter());
        PbiFilter moveAssign;
        moveAssign = tests::emptyFilter();

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 0, 4);
        tests::checkFilterInternals(moveCtor,   PbiFilter::INTERSECT, 0, 4);
        tests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 0, 4);

        EXPECT_EQ(4, original.Lookup(tests::shared_index).size());
        EXPECT_EQ(4, moveCtor.Lookup(tests::shared_index).size());
        EXPECT_EQ(4, moveAssign.Lookup(tests::shared_index).size());
    }
    { // with children
        const auto original = tests::simpleFilter();

        PbiFilter moveCtor(tests::simpleFilter());
        PbiFilter moveAssign;
        moveAssign = tests::simpleFilter();

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 1, 0);
        tests::checkFilterInternals(moveCtor,   PbiFilter::INTERSECT, 1, 0);
        tests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 1, 0);

        EXPECT_EQ(0, original.Lookup(tests::shared_index).size());
        EXPECT_EQ(0, moveCtor.Lookup(tests::shared_index).size());
        EXPECT_EQ(0, moveAssign.Lookup(tests::shared_index).size());
    }
}

TEST(PbiFilterTest, SortsAndUniquesChildFilterResultsOk)
{
    const auto childFilter = tests::SortUniqueTestFilter{ };
    const auto filter = PbiFilter{ childFilter };
    EXPECT_EQ(IndexList({ 2, 7, 0, 3, 4, 1, 8 }), childFilter.Lookup(tests::shared_index));
    EXPECT_EQ(IndexList({ 0, 1, 2, 3, 4, 7, 8 }), filter.Lookup(tests::shared_index));
}

TEST(PbiFilterTest, UnionOk)
{
    { // empty
        { // copy
            const auto emptyFilter = tests::emptyFilter();
            const auto emptyFilter2 = tests::emptyFilter();
            const auto u = PbiFilter::Union({ emptyFilter, emptyFilter2 });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, 4);
        }
        { // move
            const auto u = PbiFilter::Union({ PbiFilter{ }, PbiFilter{ } });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, 4);
        }
    }

    { // with (no-data) children - just checking composition
        { // copy
            const auto simpleFilter = tests::SimpleFilter{ };
            const auto simpleFilter2 = tests::SimpleFilter{ };
            const auto u = PbiFilter::Union({ simpleFilter, simpleFilter2 });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, 0);
        }
        { // move
            const auto u = PbiFilter::Union({ tests::SimpleFilter{ }, tests::SimpleFilter{ } });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, 0);
        }
    }

    { // 2-child union, results sorted & unique-d by PbiFilter

        const auto child1 = tests::SortUniqueTestFilter{ };
        const auto child2 = tests::SortUniqueTestFilter2{ };
        const auto u = PbiFilter::Union({ child1, child2 });
        EXPECT_EQ(IndexList({ 2, 7, 0, 3, 4, 1, 8 }),   child1.Lookup(tests::shared_index));
        EXPECT_EQ(IndexList({ 3, 7, 5 }),               child2.Lookup(tests::shared_index));
        EXPECT_EQ(IndexList({ 0, 1, 2, 3, 4, 5, 7, 8}), u.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, IntersectOk)
{
    { // empty
        { // copy
            const auto emptyFilter = tests::emptyFilter();
            const auto emptyFilter2 = tests::emptyFilter();
            const auto i = PbiFilter::Intersection({ emptyFilter, emptyFilter2 });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, 4);
        }
        { // move
            const auto i = PbiFilter::Intersection({ PbiFilter{ }, PbiFilter{ } });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, 4);
        }
    }

    { // with (no-data) children - just checking composition
        { // copy
            const auto simpleFilter = tests::SimpleFilter{ };
            const auto simpleFilter2 = tests::SimpleFilter{ };
            const auto i = PbiFilter::Intersection({ simpleFilter, simpleFilter2 });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, 0);
        }
        { // move
            const auto i = PbiFilter::Intersection({ tests::SimpleFilter{ }, tests::SimpleFilter{ } });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, 0);
        }
    }

    { // 2-child intersect, sorted & unique-d by PbiFilter

        const auto child1 = tests::SortUniqueTestFilter{ };
        const auto child2 = tests::SortUniqueTestFilter2{ };
        const auto i = PbiFilter::Intersection({ child1, child2 });
        EXPECT_EQ(IndexList({ 2, 7, 0, 3, 4, 1, 8 }), child1.Lookup(tests::shared_index));
        EXPECT_EQ(IndexList({ 3, 7, 5 }),             child2.Lookup(tests::shared_index));
        EXPECT_EQ(IndexList({ 3, 7 }),                i.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, AlignedEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4055 } };
        EXPECT_EQ(IndexList({ 1 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4055, Compare::NOT_EQUAL } };
        EXPECT_EQ(IndexList({ 0, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4000, Compare::LESS_THAN } };
        EXPECT_EQ(IndexList({ 0 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 5560, Compare::GREATER_THAN } };
        EXPECT_EQ(IndexList({ 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 5560, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 2, 3 }), filter.Lookup(tests::shared_index));
    }

    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 7000, Compare::GREATER_THAN } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, AlignedLengthFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedLengthFilter{ 500, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedLengthFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 2 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, AlignedStartFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 2600, Compare::LESS_THAN } };
        EXPECT_EQ(IndexList({ 0, 1 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 4102, Compare::GREATER_THAN } };
        EXPECT_EQ(IndexList({ 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 4102, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 6000, Compare::GREATER_THAN } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, AlignedStrandFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::FORWARD } };
        EXPECT_EQ(IndexList({ 0, 2 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::REVERSE } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::FORWARD, Compare::NOT_EQUAL } }; // same as Strand::REVERSE
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }

    // unsupported compare types throw
    EXPECT_THROW(PbiAlignedStrandFilter(Strand::FORWARD, Compare::LESS_THAN),          std::runtime_error);
    EXPECT_THROW(PbiAlignedStrandFilter(Strand::FORWARD, Compare::LESS_THAN_EQUAL),    std::runtime_error);
    EXPECT_THROW(PbiAlignedStrandFilter(Strand::FORWARD, Compare::GREATER_THAN),       std::runtime_error);
    EXPECT_THROW(PbiAlignedStrandFilter(Strand::FORWARD, Compare::GREATER_THAN_EQUAL), std::runtime_error);
}

TEST(PbiFilterTest, BarcodeFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeFilter{ 17 } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeFilter{ 18 } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeFilter{ 0 } };
        EXPECT_EQ(IndexList({ 0 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, BarcodeForwardFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ 17 } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ 400 } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ {0, 256} } };
        EXPECT_EQ(IndexList({ 0, 2 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, BarcodeQualityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeQualityFilter{ 80, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeQualityFilter{ 40, Compare::LESS_THAN } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, BarcodeReverseFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ 18 } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ 400 } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ {1, 257} } };
        EXPECT_EQ(IndexList({ 0, 2 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, BarcodesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ 17, 18 } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ 17, 19 } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ std::make_pair(17,18) } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, IdentityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiIdentityFilter{ 0.95, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 3 }), filter.Lookup(tests::shared_index));
    }
}

//TEST(PbiFilterTest, LocalContextFlagFilterOk)
//{
//    auto f = PbiLocalContextFlagFilter{ LocalContextFlags::ADAPTER_AFTER };
//    (void)f;
//    EXPECT_TRUE(false);
//}

TEST(PbiFilterTest, MapQualityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiMapQualityFilter{ 254 } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiMapQualityFilter{ 254, Compare::NOT_EQUAL } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, MovieNameFilterOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/test_group_query/test2.bam" } };
    const auto index = PbiIndex{ bamFile.PacBioIndexFilename() };

    {
        const auto filter = PbiFilter{ PbiMovieNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0" } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(index));
    }
    {
        const auto filter = PbiFilter{ PbiMovieNameFilter{ "does_not_exist" } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(index));
    }
    {
        const auto names = vector<string>{"does_not_exist",
                                          "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0"};
        const auto filter = PbiFilter{ PbiMovieNameFilter{ names } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(index));
    }
}

TEST(PbiFilterTest, NumDeletedBasesFilterOk)
{
    // del: { 12, 38, 45, 11} - calculated from raw data, not stored directly in testing object or read from PBI file

    {
        const auto filter = PbiFilter{ PbiNumDeletedBasesFilter{ 12, Compare::LESS_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 0, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiNumDeletedBasesFilter{ 45, Compare::EQUAL } };
        EXPECT_EQ(IndexList({ 2 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, NumInsertedBasesFilterOk)
{
    // ins: { 17, 63, 65, 20 }  - calculated from raw data, not stored directly testing object or read from PBI file

    {
        const auto filter = PbiFilter{ PbiNumInsertedBasesFilter{ 63, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 2 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiNumInsertedBasesFilter{ 17, Compare::NOT_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, NumMatchesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiNumMatchesFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 2 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiNumMatchesFilter{ 400, Compare::LESS_THAN } };
        EXPECT_EQ(IndexList({ 0 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, NumMismatchesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiNumMismatchesFilter{ 0, Compare::EQUAL } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiNumMismatchesFilter{ 0, Compare::NOT_EQUAL } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, QueryEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryEndFilter{ 4055 } };
        EXPECT_EQ(IndexList({ 1 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiQueryEndFilter{ 6200, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, QueryLengthFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryLengthFilter{ 500, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiQueryLengthFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 1, 2 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, QueryNameFilterOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/test_group_query/test2.bam" } };
    const auto index = PbiIndex{ bamFile.PacBioIndexFilename() };

    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055" } };
        EXPECT_EQ(IndexList({ 1 }), filter.Lookup(index));
    }
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237" } };
        EXPECT_EQ(IndexList({ 3 }), filter.Lookup(index));
    }

    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "does_not_exist/0/0_0" } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(index));
    }
    {
        const auto names = vector<string>{"m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055",
                                          "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"};
        const auto filter = PbiFilter{ PbiQueryNameFilter{ names } };
        EXPECT_EQ(IndexList({ 1, 3 }), filter.Lookup(index));
    }

    // invalid QNAME syntax throws
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "" } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(index));
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo" } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(index));
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo/bar" } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(index));
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo/bar/baz_bam" } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(index));
    },
    std::runtime_error);
}

TEST(PbiFilterTest, QueryStartFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 4101 } };
        EXPECT_EQ(IndexList({ 2 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 5000 } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 5000, Compare::GREATER_THAN } };
        EXPECT_EQ(IndexList({ 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, ReadAccuracyFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReadAccuracyFilter{ 0.9 } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiReadAccuracyFilter{ 0.9, Compare::GREATER_THAN } };
        EXPECT_EQ(IndexList({ 0, 2 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, ReadGroupFilterOk)
{
    { // numeric ID
        const auto filter = PbiReadGroupFilter{ -1197849594 };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));

        const auto filter2 = PbiReadGroupFilter{ 200 };
        EXPECT_EQ(IndexList({ }), filter2.Lookup(tests::shared_index));
    }
    { // string ID
        const auto filter = PbiReadGroupFilter{ "b89a4406" };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));

        const auto filter2 = PbiReadGroupFilter{ "b89a4406" };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter2.Lookup(tests::shared_index));
    }
    { // ReadGroupInfo object
        const auto rg = ReadGroupInfo{ "b89a4406" };
        const auto filter = PbiReadGroupFilter{ rg };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    { // multi-ID
        const auto ids = vector<int32_t>({-1197849594, 200});
        const auto filter = PbiReadGroupFilter{ ids };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    { // multi-string
        const auto ids = vector<string>({"b89a4406", "deadbeef"});
        const auto filter = PbiReadGroupFilter{ ids };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    { // multi-ReadGroupInfo
        const auto ids = vector<ReadGroupInfo>({ ReadGroupInfo("b89a4406"), ReadGroupInfo("deadbeef")});
        const auto filter = PbiReadGroupFilter{ ids };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, ReferenceEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReferenceEndFilter{ 9900 } };
        EXPECT_EQ(IndexList({ 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiReferenceEndFilter{ 9900, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 0, 1, 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, ReferenceIdFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ 0 } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ 0, Compare::NOT_EQUAL } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
    {
        const auto ids = vector<int32_t>({0, 42});
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ ids } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, ReferenceNameFilterOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/test_group_query/test2.bam" } };
    const auto index = PbiIndex{ bamFile.PacBioIndexFilename() };

    {
        const auto filter = PbiFilter{ PbiReferenceNameFilter{ "lambda_NEB3011" } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(index));
    }
    {
        const auto filter = PbiFilter{ PbiReferenceNameFilter{ "lambda_NEB3011", Compare::NOT_EQUAL } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(index));
    }
    {
        const auto names = vector<string>({ "lambda_NEB3011" }); // this file only has 1 :(
        const auto filter = PbiFilter{ PbiReferenceNameFilter{ names } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(index));
    }

    // unsupported compare types throw
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::LESS_THAN),          std::runtime_error);
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::LESS_THAN_EQUAL),    std::runtime_error);
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::GREATER_THAN),       std::runtime_error);
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::GREATER_THAN_EQUAL), std::runtime_error);
}

TEST(PbiFilterTest, ReferenceStartFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReferenceStartFilter{ 8453 } };
        EXPECT_EQ(IndexList({ 1 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiReferenceStartFilter{ 9200, Compare::GREATER_THAN_EQUAL } };
        EXPECT_EQ(IndexList({ 0, 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, ZmwFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiZmwFilter{ 14743 } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
    {
        const auto filter = PbiFilter{ PbiZmwFilter{ 14743, Compare::NOT_EQUAL } };
        EXPECT_EQ(IndexList({ }), filter.Lookup(tests::shared_index));
    }
    {
        const auto zmws = vector<int32_t>({14743,42,200});
        const auto filter = PbiFilter{ PbiZmwFilter{ zmws } };
        EXPECT_EQ(IndexList({ 0, 1, 2, 3 }), filter.Lookup(tests::shared_index));
    }
}

TEST(PbiFilterTest, FromDataSetOk)
{
    const auto expectedFilter =
        PbiFilter::Union(
        {
            PbiFilter::Intersection(
            {
                PbiZmwFilter{ 14743 },
                PbiReadAccuracyFilter { 0.9, Compare::GREATER_THAN_EQUAL }
            }),

            PbiReferenceStartFilter { 9200, Compare::GREATER_THAN_EQUAL }
        });


    auto properties1 = Properties{ };
    properties1.Add(Property{ "zm", "14743",  "==" });
    properties1.Add(Property{ "rq", "0.9", ">=" });

    auto datasetFilter1 = Filter{ };
    datasetFilter1.Properties(properties1);

    auto properties2 = Properties{ };
    properties2.Add(Property{ "pos", "9200", ">=" });

    auto datasetFilter2 = Filter{ };
    datasetFilter2.Properties(properties2);

    auto datasetFilters = Filters{ };
    datasetFilters.Add(datasetFilter1);
    datasetFilters.Add(datasetFilter2);
    auto dataset = DataSet{ };
    dataset.Filters(datasetFilters);

    const auto generatedFilter = PbiFilter::FromDataSet(dataset);
    EXPECT_EQ(expectedFilter.Lookup(tests::shared_index),
              generatedFilter.Lookup(tests::shared_index));
}
