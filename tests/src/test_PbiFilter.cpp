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

static const PbiRawData shared_index = test2Bam_RawIndex();

static
void checkFilterRows(const PbiFilter& filter, const std::vector<size_t> expectedRows)
{
    for (size_t row : expectedRows)
        EXPECT_TRUE(filter.Accepts(shared_index, row));
}

static
void checkFilterInternals(const PbiFilter& filter,
                          const PbiFilter::CompositionType expectedType,
                          const size_t expectedNumChildren,
                          const std::vector<size_t> expectedRows)
{
    EXPECT_EQ(expectedType,        filter.d_->type_);
    EXPECT_EQ(expectedNumChildren, filter.d_->filters_.size());
    checkFilterRows(filter, expectedRows);
}

struct SimpleFilter
{
    bool Accepts(const PbiRawData& idx, const size_t row) const
    { (void)idx; (void)row; return true; }
};

struct NoncompliantFilter { };

struct SortUniqueTestFilter
{
    bool Accepts(const PbiRawData& idx, const size_t row) const
    {
        (void)idx;
        switch(row) {
            case 0: // fall through
            case 1: // .
            case 2: // .
            case 3: // .
            case 4: // .
            case 7: // .
            case 8: return true;
            default:
                return false;
        }
    }
};

struct SortUniqueTestFilter2
{
    bool Accepts(const PbiRawData& idx, const size_t row) const
    {
        (void)idx;
        switch(row) {
            case 3: // fall through
            case 7: // .
            case 5: return true;
            default:
                return false;
        }
    }
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
    tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
}

TEST(PbiFilterTest, CompositionOk)
{
    auto filter = PbiFilter{ };
    filter.Add(PbiFilter{ });
    tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
}

TEST(PbiFilterTest, CustomFilterOk)
{
    { // ctor
        auto filter = PbiFilter{ tests::SimpleFilter{ } };
        tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, std::vector<size_t>{});
    }
    { // Add
        auto filter = PbiFilter{ };
        filter.Add(tests::SimpleFilter{ });
        tests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, std::vector<size_t>{});
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

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        tests::checkFilterInternals(copyCtor,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        tests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
    }
    { // with children
        const auto original = PbiFilter{ tests::SimpleFilter{ } };

        PbiFilter copyCtor(original);
        PbiFilter copyAssign;
        copyAssign = original;

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 1, std::vector<size_t>{});
        tests::checkFilterInternals(copyCtor,   PbiFilter::INTERSECT, 1, std::vector<size_t>{});
        tests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 1, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, MoveOk)
{
    { // empty
        const auto original = tests::emptyFilter();

        PbiFilter moveCtor(tests::emptyFilter());
        PbiFilter moveAssign;
        moveAssign = tests::emptyFilter();

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        tests::checkFilterInternals(moveCtor,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        tests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
    }
    { // with children
        const auto original = tests::simpleFilter();

        PbiFilter moveCtor(tests::simpleFilter());
        PbiFilter moveAssign;
        moveAssign = tests::simpleFilter();

        tests::checkFilterInternals(original,   PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
        tests::checkFilterInternals(moveCtor,   PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
        tests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
    }
}

TEST(PbiFilterTest, SortsAndUniquesChildFilterResultsOk)
{
    const auto childFilter = tests::SortUniqueTestFilter{ };
    const auto filter = PbiFilter{ childFilter };
    tests::checkFilterRows(childFilter, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
    tests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3, 4, 7, 8});
}

TEST(PbiFilterTest, UnionOk)
{
    { // empty
        { // copy
            const auto emptyFilter = tests::emptyFilter();
            const auto emptyFilter2 = tests::emptyFilter();
            const auto u = PbiFilter::Union({ emptyFilter, emptyFilter2 });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0,1,2,3});
        }
        { // move
            const auto u = PbiFilter::Union({ PbiFilter{ }, PbiFilter{ } });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0,1,2,3});
        }
    }

    { // with (no-data) children - just checking composition
        { // copy
            const auto simpleFilter = tests::SimpleFilter{ };
            const auto simpleFilter2 = tests::SimpleFilter{ };
            const auto u = PbiFilter::Union({ simpleFilter, simpleFilter2 });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{});
        }
        { // move
            const auto u = PbiFilter::Union({ tests::SimpleFilter{ }, tests::SimpleFilter{ } });
            tests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{});
        }
    }

    { // 2-child union, results sorted & unique-d by PbiFilter

        const auto child1 = tests::SortUniqueTestFilter{ };
        const auto child2 = tests::SortUniqueTestFilter2{ };
        const auto u = PbiFilter::Union({ child1, child2 });

        tests::checkFilterRows(child1, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
        tests::checkFilterRows(child2, std::vector<size_t>{3, 7, 5});
        tests::checkFilterRows(u, std::vector<size_t>{0, 1, 2, 3, 4, 5, 7, 8});
    }
}

TEST(PbiFilterTest, IntersectOk)
{
    { // empty
        { // copy
            const auto emptyFilter = tests::emptyFilter();
            const auto emptyFilter2 = tests::emptyFilter();
            const auto i = PbiFilter::Intersection({ emptyFilter, emptyFilter2 });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{0,1,2,3});
        }
        { // move
            const auto i = PbiFilter::Intersection({ PbiFilter{ }, PbiFilter{ } });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{0,1,2,3});
        }
    }

    { // with (no-data) children - just checking composition
        { // copy
            const auto simpleFilter = tests::SimpleFilter{ };
            const auto simpleFilter2 = tests::SimpleFilter{ };
            const auto i = PbiFilter::Intersection({ simpleFilter, simpleFilter2 });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{});
        }
        { // move
            const auto i = PbiFilter::Intersection({ tests::SimpleFilter{ }, tests::SimpleFilter{ } });
            tests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{});
        }
    }

    { // 2-child intersect, sorted & unique-d by PbiFilter

        const auto child1 = tests::SortUniqueTestFilter{ };
        const auto child2 = tests::SortUniqueTestFilter2{ };
        const auto i = PbiFilter::Intersection({ child1, child2 });

        tests::checkFilterRows(child1, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
        tests::checkFilterRows(child2, std::vector<size_t>{3, 7, 5 });
        tests::checkFilterRows(i, std::vector<size_t>{3, 7});
    }
}

TEST(PbiFilterTest, AlignedEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4055 } };
        tests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4055, Compare::NOT_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4000, Compare::LESS_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{0});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 5560, Compare::GREATER_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 5560, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{2,3});
    }

    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 7000, Compare::GREATER_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, AlignedLengthFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedLengthFilter{ 500, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedLengthFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
}

TEST(PbiFilterTest, AlignedStartFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 2600, Compare::LESS_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 4102, Compare::GREATER_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 4102, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{2,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 6000, Compare::GREATER_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{ });
    }
}

TEST(PbiFilterTest, AlignedStrandFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::FORWARD } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::REVERSE } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::FORWARD, Compare::NOT_EQUAL } }; // same as Strand::REVERSE
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
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
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeFilter{ 18 } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeFilter{ 0 } };
        tests::checkFilterRows(filter, std::vector<size_t>{0});
    }
}

TEST(PbiFilterTest, BarcodeForwardFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ 17 } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ 400 } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ {0, 256} } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
}

TEST(PbiFilterTest, BarcodeQualityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeQualityFilter{ 80, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeQualityFilter{ 40, Compare::LESS_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, BarcodeReverseFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ 18 } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ 400 } };
        tests::checkFilterRows(filter, std::vector<size_t>{ });
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ {1, 257} } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
}

TEST(PbiFilterTest, BarcodesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ 17, 18 } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ 17, 19 } };
        tests::checkFilterRows(filter, std::vector<size_t>{ });
    }
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ std::make_pair(17,18) } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
}

TEST(PbiFilterTest, IdentityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiIdentityFilter{ 0.95, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{3});
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
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiMapQualityFilter{ 254, Compare::NOT_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, MovieNameFilterOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/test_group_query/test2.bam" } };
    const auto index = PbiRawData{ bamFile.PacBioIndexFilename() };

    {
        const auto filter = PbiFilter{ PbiMovieNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0" } };
        const auto expectedRows = std::vector<size_t>{0,1,2,3};
        for (size_t row : expectedRows)
            EXPECT_TRUE(filter.Accepts(index, row));
    }
    {
        const auto filter = PbiFilter{ PbiMovieNameFilter{ "does_not_exist" } };
        const auto expectedRows = std::vector<size_t>{};
        for (size_t row : expectedRows)
            EXPECT_TRUE(filter.Accepts(index, row));
    }
    {
        const auto names = vector<string>{"does_not_exist",
                                          "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0"};
        const auto filter = PbiFilter{ PbiMovieNameFilter{ names } };
        const auto expectedRows = std::vector<size_t>{0,1,2,3};
        for (size_t row : expectedRows)
            EXPECT_TRUE(filter.Accepts(index, row));
    }
}

TEST(PbiFilterTest, NumDeletedBasesFilterOk)
{
    // del: { 12, 38, 45, 11} - calculated from raw data, not stored directly in testing object or read from PBI file

    {
        const auto filter = PbiFilter{ PbiNumDeletedBasesFilter{ 12, Compare::LESS_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,3});
    }
    {
        const auto filter = PbiFilter{ PbiNumDeletedBasesFilter{ 45, Compare::EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{2});
    }
}

TEST(PbiFilterTest, NumInsertedBasesFilterOk)
{
    // ins: { 17, 63, 65, 20 }  - calculated from raw data, not stored directly testing object or read from PBI file

    {
        const auto filter = PbiFilter{ PbiNumInsertedBasesFilter{ 63, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
    {
        const auto filter = PbiFilter{ PbiNumInsertedBasesFilter{ 17, Compare::NOT_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,2,3});
    }
}

TEST(PbiFilterTest, NumMatchesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiNumMatchesFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
    {
        const auto filter = PbiFilter{ PbiNumMatchesFilter{ 400, Compare::LESS_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{0});
    }
}

TEST(PbiFilterTest, NumMismatchesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiNumMismatchesFilter{ 0, Compare::EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiNumMismatchesFilter{ 0, Compare::NOT_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, QueryEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryEndFilter{ 4055 } };
        tests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiQueryEndFilter{ 6200, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(PbiFilterTest, QueryLengthFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryLengthFilter{ 500, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiQueryLengthFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
}

TEST(PbiFilterTest, QueryNameFilterOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/test_group_query/test2.bam" } };
    const auto index = PbiIndex{ bamFile.PacBioIndexFilename() };

    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055" } };
        tests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237" } };
        tests::checkFilterRows(filter, std::vector<size_t>{3});
    }

    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "does_not_exist/0/0_0" } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto names = vector<string>{"m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055",
                                          "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"};
        const auto filter = PbiFilter{ PbiQueryNameFilter{ names } };
        tests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }

    // invalid QNAME syntax throws
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "" } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo" } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo/bar" } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo/bar/baz_bam" } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::exception); // come back to see why this is not runtime_error but something else
}

TEST(PbiFilterTest, QueryStartFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 4101 } };
        tests::checkFilterRows(filter, std::vector<size_t>{2});
    }
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 5000 } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 5000, Compare::GREATER_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(PbiFilterTest, ReadAccuracyFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReadAccuracyFilter{ 0.9 } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{ PbiReadAccuracyFilter{ 0.9, Compare::GREATER_THAN } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
}

TEST(PbiFilterTest, ReadGroupFilterOk)
{
    { // numeric ID
        const auto filter = PbiReadGroupFilter{ -1197849594 };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});

        const auto filter2 = PbiReadGroupFilter{ 200 };
        tests::checkFilterRows(filter2, std::vector<size_t>{});
    }
    { // string ID
        const auto filter = PbiReadGroupFilter{ "b89a4406" };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});

        const auto filter2 = PbiReadGroupFilter{ "b89a4406" };
        tests::checkFilterRows(filter2, std::vector<size_t>{0,1,2,3});
    }
    { // ReadGroupInfo object
        const auto rg = ReadGroupInfo{ "b89a4406" };
        const auto filter = PbiReadGroupFilter{ rg };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    { // multi-ID
        const auto ids = vector<int32_t>({-1197849594, 200});
        const auto filter = PbiReadGroupFilter{ ids };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    { // multi-string
        const auto ids = vector<string>({"b89a4406", "deadbeef"});
        const auto filter = PbiReadGroupFilter{ ids };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    { // multi-ReadGroupInfo
        const auto ids = vector<ReadGroupInfo>({ ReadGroupInfo("b89a4406"), ReadGroupInfo("deadbeef")});
        const auto filter = PbiReadGroupFilter{ ids };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
}

TEST(PbiFilterTest, ReferenceEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReferenceEndFilter{ 9900 } };
        tests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{ PbiReferenceEndFilter{ 9900, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,3});
    }
}

TEST(PbiFilterTest, ReferenceIdFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ 0 } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ 0, Compare::NOT_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto ids = vector<int32_t>({0, 42});
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ ids } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
}

TEST(PbiFilterTest, ReferenceNameFilterOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/test_group_query/test2.bam" } };
    const auto index = PbiRawData{ bamFile.PacBioIndexFilename() };

    {
        const auto filter = PbiFilter{ PbiReferenceNameFilter{ "lambda_NEB3011" } };
        const auto expectedRows = std::vector<size_t>{0,1,2,3};
        for (size_t row : expectedRows)
            EXPECT_TRUE(filter.Accepts(index, row));

    }
    {
        const auto filter = PbiFilter{ PbiReferenceNameFilter{ "lambda_NEB3011", Compare::NOT_EQUAL } };
        const auto expectedRows = std::vector<size_t>{};
        for (size_t row : expectedRows)
            EXPECT_TRUE(filter.Accepts(index, row));
    }
    {
        const auto names = vector<string>({ "lambda_NEB3011" }); // this file only has 1 :(
        const auto filter = PbiFilter{ PbiReferenceNameFilter{ names } };
        const auto expectedRows = std::vector<size_t>{0,1,2,3};
        for (size_t row : expectedRows)
            EXPECT_TRUE(filter.Accepts(index, row));
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
        tests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiReferenceStartFilter{ 9200, Compare::GREATER_THAN_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,3});
    }
}

TEST(PbiFilterTest, ZmwFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiZmwFilter{ 14743 } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiZmwFilter{ 14743, Compare::NOT_EQUAL } };
        tests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto zmws = vector<int32_t>({14743,42,200});
        const auto filter = PbiFilter{ PbiZmwFilter{ zmws } };
        tests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
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

    for (size_t i = 0; i < tests::shared_index.NumReads(); ++i) {
        EXPECT_EQ(expectedFilter.Accepts(tests::shared_index, i),
                  generatedFilter.Accepts(tests::shared_index, i));
    }
}
