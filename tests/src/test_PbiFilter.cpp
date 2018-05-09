// Author: Derek Barnett

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/PbiFilter.h>

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;

namespace PbiFilterTests {

// helper structs & methods

static
PbiRawData test2Bam_RawIndex()
{
    PbiRawData index;
    index.NumReads(4);

    PbiRawBasicData& subreadData = index.BasicData();
    subreadData.rgId_       = { -1197849594, -1197849594, -1197849594, -1197849594 };
    subreadData.qStart_     = { 2114, 2579, 4101, 5615 };
    subreadData.qEnd_       = { 2531, 4055, 5571, 6237 };
    subreadData.holeNumber_ = { 14743, 14743, 14743, 14743 };
    subreadData.readQual_   = { 0.901, 0.601, 0.901, 0.601 };
    subreadData.ctxtFlag_   = { 0, 1, 2, 3 };
    subreadData.fileOffset_ = { 35651584, 35655125, 35667128, 35679170 };

    PbiRawMappedData& mappedData = index.MappedData();
    mappedData.tId_       = { 0, 0, 0, 0 };
    mappedData.tStart_    = { 9507, 8453, 8455, 9291 };
    mappedData.tEnd_      = { 9903, 9902, 9893, 9900 };
    mappedData.aStart_    = { 2130, 2581, 4102, 5619 };
    mappedData.aEnd_      = { 2531, 4055, 5560, 6237 };
    mappedData.revStrand_ = { 0, 1, 0, 1 };
    mappedData.mapQV_     = { 254, 254, 254, 254 };
    mappedData.nM_        = { 384, 1411, 1393, 598 };
    mappedData.nMM_       = { 0, 0, 0, 0 };

    PbiRawBarcodeData& barcodeData = index.BarcodeData();
    barcodeData.bcForward_ = { 0, 17, 256, 17 };
    barcodeData.bcReverse_ = { 1, 18, 257, 18 };
    barcodeData.bcQual_    = { 42, 80, 42, 110 };

    PbiRawReferenceData& referenceData = index.ReferenceData();
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
    EXPECT_EQ(expectedType,        filter.Type());
    EXPECT_EQ(expectedNumChildren, filter.NumChildren());
    checkFilterRows(filter, expectedRows);
}

struct SimpleFilter
{
    bool Accepts(const PbiRawData& /* idx */, const size_t /* row */) const
    { /*()idx; ()row;*/ return true; }
};

struct NoncompliantFilter { };

struct SortUniqueTestFilter
{
    bool Accepts(const PbiRawData& /* idx */, const size_t row) const
    {
//        ()idx;
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
    bool Accepts(const PbiRawData& /* idx */, const size_t row) const
    {
//        ()idx;
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
PbiFilter emptyFilter()
{ return PbiFilter{ }; }

static inline
PbiFilter simpleFilter()
{ return PbiFilter{ SimpleFilter{ } }; }

} // namespace PbiFilterTests

TEST(PbiFilterTest, DefaultCtorOk)
{
    auto filter = PbiFilter{ };
    PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
}

TEST(PbiFilterTest, CompositionOk)
{
    auto filter = PbiFilter{ };
    filter.Add(PbiFilter{ });
    PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
}

TEST(PbiFilterTest, CustomFilterOk)
{
    { // ctor
        auto filter = PbiFilter{ PbiFilterTests::SimpleFilter{ } };
        PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, std::vector<size_t>{});
    }
    { // Add
        auto filter = PbiFilter{ };
        filter.Add(PbiFilterTests::SimpleFilter{ });
        PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1, std::vector<size_t>{});
    }

//    PbiFilter shouldNotCompile = PbiFilter{ PbiFilterTests::NoncompliantFilter{ } };                       // <-- when uncommented, should not compile
//    PbiFilter shouldNotCompileEither; shouldNotCompileEither.Add(PbiFilterTests::NoncompliantFilter{ });   // <-- when uncommented, should not compile
}

TEST(PbiFilterTest, CopyOk)
{
    { // empty
        const auto original = PbiFilter{ };

        PbiFilter copyCtor(original);
        PbiFilter copyAssign;
        copyAssign = original;

        PbiFilterTests::checkFilterInternals(original,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        PbiFilterTests::checkFilterInternals(copyCtor,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        PbiFilterTests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
    }
    { // with children
        const auto original = PbiFilter{ PbiFilterTests::SimpleFilter{ } };

        PbiFilter copyCtor(original);
        PbiFilter copyAssign;
        copyAssign = original;

        PbiFilterTests::checkFilterInternals(original,   PbiFilter::INTERSECT, 1, std::vector<size_t>{});
        PbiFilterTests::checkFilterInternals(copyCtor,   PbiFilter::INTERSECT, 1, std::vector<size_t>{});
        PbiFilterTests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 1, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, MoveOk)
{
    { // empty
        const auto original = PbiFilterTests::emptyFilter();

        PbiFilter moveCtor(PbiFilterTests::emptyFilter());
        PbiFilter moveAssign;
        moveAssign = PbiFilterTests::emptyFilter();

        PbiFilterTests::checkFilterInternals(original,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        PbiFilterTests::checkFilterInternals(moveCtor,   PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
        PbiFilterTests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 0, std::vector<size_t>{0,1,2,3});
    }
    { // with children
        const auto original = PbiFilterTests::simpleFilter();

        PbiFilter moveCtor(PbiFilterTests::simpleFilter());
        PbiFilter moveAssign;
        moveAssign = PbiFilterTests::simpleFilter();

        PbiFilterTests::checkFilterInternals(original,   PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
        PbiFilterTests::checkFilterInternals(moveCtor,   PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
        PbiFilterTests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 1, std::vector<size_t>{0,1,2,3});
    }
}

TEST(PbiFilterTest, SortsAndUniquesChildFilterResultsOk)
{
    const auto childFilter = PbiFilterTests::SortUniqueTestFilter{ };
    const auto filter = PbiFilter{ childFilter };
    PbiFilterTests::checkFilterRows(childFilter, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
    PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3, 4, 7, 8});
}

TEST(PbiFilterTest, UnionOk)
{
    { // empty
        { // copy
            const auto emptyFilter = PbiFilterTests::emptyFilter();
            const auto emptyFilter2 = PbiFilterTests::emptyFilter();
            const auto u = PbiFilter::Union({ emptyFilter, emptyFilter2 });
            PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0,1,2,3});
        }
        { // move
            const auto u = PbiFilter::Union({ PbiFilter{ }, PbiFilter{ } });
            PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0,1,2,3});
        }
    }

    { // with (no-data) children - just checking composition
        { // copy
            const auto simpleFilter = PbiFilterTests::SimpleFilter{ };
            const auto simpleFilter2 = PbiFilterTests::SimpleFilter{ };
            const auto u = PbiFilter::Union({ simpleFilter, simpleFilter2 });
            PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{});
        }
        { // move
            const auto u = PbiFilter::Union({ PbiFilterTests::SimpleFilter{ }, PbiFilterTests::SimpleFilter{ } });
            PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{});
        }
    }

    { // 2-child union, results sorted & unique-d by PbiFilter

        const auto child1 = PbiFilterTests::SortUniqueTestFilter{ };
        const auto child2 = PbiFilterTests::SortUniqueTestFilter2{ };
        const auto u = PbiFilter::Union({ child1, child2 });

        PbiFilterTests::checkFilterRows(child1, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
        PbiFilterTests::checkFilterRows(child2, std::vector<size_t>{3, 7, 5});
        PbiFilterTests::checkFilterRows(u, std::vector<size_t>{0, 1, 2, 3, 4, 5, 7, 8});
    }
}

TEST(PbiFilterTest, IntersectOk)
{
    { // empty
        { // copy
            const auto emptyFilter = PbiFilterTests::emptyFilter();
            const auto emptyFilter2 = PbiFilterTests::emptyFilter();
            const auto i = PbiFilter::Intersection({ emptyFilter, emptyFilter2 });
            PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{0,1,2,3});
        }
        { // move
            const auto i = PbiFilter::Intersection({ PbiFilter{ }, PbiFilter{ } });
            PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{0,1,2,3});
        }
    }

    { // with (no-data) children - just checking composition
        { // copy
            const auto simpleFilter = PbiFilterTests::SimpleFilter{ };
            const auto simpleFilter2 = PbiFilterTests::SimpleFilter{ };
            const auto i = PbiFilter::Intersection({ simpleFilter, simpleFilter2 });
            PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{});
        }
        { // move
            const auto i = PbiFilter::Intersection({ PbiFilterTests::SimpleFilter{ }, PbiFilterTests::SimpleFilter{ } });
            PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{});
        }
    }

    { // 2-child intersect, sorted & unique-d by PbiFilter

        const auto child1 = PbiFilterTests::SortUniqueTestFilter{ };
        const auto child2 = PbiFilterTests::SortUniqueTestFilter2{ };
        const auto i = PbiFilter::Intersection({ child1, child2 });

        PbiFilterTests::checkFilterRows(child1, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
        PbiFilterTests::checkFilterRows(child2, std::vector<size_t>{3, 7, 5 });
        PbiFilterTests::checkFilterRows(i, std::vector<size_t>{3, 7});
    }
}

TEST(PbiFilterTest, AlignedEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4055 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4055, Compare::NOT_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 4000, Compare::LESS_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 5560, Compare::GREATER_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 5560, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2,3});
    }

    {
        const auto filter = PbiFilter{ PbiAlignedEndFilter{ 7000, Compare::GREATER_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, AlignedLengthFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedLengthFilter{ 500, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedLengthFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
}

TEST(PbiFilterTest, AlignedStartFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 2600, Compare::LESS_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 4102, Compare::GREATER_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 4102, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStartFilter{ 6000, Compare::GREATER_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{ });
    }
}

TEST(PbiFilterTest, AlignedStrandFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::FORWARD } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::REVERSE } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiAlignedStrandFilter{ Strand::FORWARD, Compare::NOT_EQUAL } }; // same as Strand::REVERSE
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
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
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeFilter{ 18 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeFilter{ 0 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
}

TEST(PbiFilterTest, BarcodeForwardFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ 17 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ 400 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeForwardFilter{ {0, 256} } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
}

TEST(PbiFilterTest, BarcodeQualityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeQualityFilter{ 80, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeQualityFilter{ 40, Compare::LESS_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, BarcodeReverseFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ 18 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ 400 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{ });
    }
    {
        const auto filter = PbiFilter{ PbiBarcodeReverseFilter{ {1, 257} } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
}

TEST(PbiFilterTest, BarcodesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ 17, 18 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ 17, 19 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{ });
    }
    {
        const auto filter = PbiFilter{ PbiBarcodesFilter{ std::make_pair(17,18) } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
}

TEST(PbiFilterTest, IdentityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiIdentityFilter{ 0.95, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(PbiFilterTest, LocalContextFilterOk)
{
    { // == NO_LOCAL_CONTEXT
        const auto filter = PbiFilter { PbiLocalContextFilter{ LocalContextFlags::NO_LOCAL_CONTEXT } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
    { // != ADAPTER_BEFORE (exact match)
        const auto filter = PbiFilter { PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,2,3});
    }
    { // contains ADAPTER_BEFORE
        const auto filter = PbiFilter { PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }
    { // does not contain ADAPTER_BEFORE
        const auto filter = PbiFilter { PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
    { // include both ADAPTER_BEFORE and ADAPTER_AFTER
        const auto filter = PbiFilter::Intersection(
        {
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS },
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::CONTAINS }
        });
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    { // exclude both ADAPTER_BEFORE and ADAPTER_AFTER
        const auto filter = PbiFilter::Intersection(
        {
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS },
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::NOT_CONTAINS }
        });
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
    { // include everything with either ADAPTER_BEFORE or ADAPTER_AFTER
        const auto filter = PbiFilter::Union(
        {
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS },
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::CONTAINS }
        });
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2,3});
    }
    { // include everything with either ADAPTER_BEFORE or ADAPTER_AFTER, but not both
        const auto filter = PbiFilter::Intersection(
        {
                PbiLocalContextFilter{ LocalContextFlags::NO_LOCAL_CONTEXT, Compare::NOT_EQUAL },
                PbiFilter::Union(
                {
                    PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS },
                    PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::NOT_CONTAINS }
                })
        });
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
}

TEST(PbiFilterTest, MapQualityFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiMapQualityFilter{ 254 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiMapQualityFilter{ 254, Compare::NOT_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, MovieNameFilterOk)
{
    const auto bamFile = BamFile{ PbbamTestsConfig::Data_Dir + std::string{ "/group/test2.bam" } };
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
        const auto names = std::vector<std::string>{"does_not_exist",
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
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,3});
    }
    {
        const auto filter = PbiFilter{ PbiNumDeletedBasesFilter{ 45, Compare::EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2});
    }
}

TEST(PbiFilterTest, NumInsertedBasesFilterOk)
{
    // ins: { 17, 63, 65, 20 }  - calculated from raw data, not stored directly testing object or read from PBI file

    {
        const auto filter = PbiFilter{ PbiNumInsertedBasesFilter{ 63, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
    {
        const auto filter = PbiFilter{ PbiNumInsertedBasesFilter{ 17, Compare::NOT_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2,3});
    }
}

TEST(PbiFilterTest, NumMatchesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiNumMatchesFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
    {
        const auto filter = PbiFilter{ PbiNumMatchesFilter{ 400, Compare::LESS_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
}

TEST(PbiFilterTest, NumMismatchesFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiNumMismatchesFilter{ 0, Compare::EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiNumMismatchesFilter{ 0, Compare::NOT_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(PbiFilterTest, QueryEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryEndFilter{ 4055 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiQueryEndFilter{ 6200, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(PbiFilterTest, QueryLengthFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryLengthFilter{ 500, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiQueryLengthFilter{ 1000, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,2});
    }
}

TEST(PbiFilterTest, QueryNameFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055" } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237" } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }

    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "does_not_exist/0/0_0" } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto names = std::vector<std::string>{"m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055",
                                          "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"};
        const auto filter = PbiFilter{ PbiQueryNameFilter{ names } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1,3});
    }

    // invalid QNAME syntax throws
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "" } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo" } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo/bar" } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::runtime_error);
    EXPECT_THROW(
    {
        const auto filter = PbiFilter{ PbiQueryNameFilter{ "foo/bar/baz_bam" } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    },
    std::exception); // come back to see why this is not runtime_error but something else
}

TEST(PbiFilterTest, QueryStartFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 4101 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2});
    }
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 5000 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{ PbiQueryStartFilter{ 5000, Compare::GREATER_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(PbiFilterTest, ReadAccuracyFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReadAccuracyFilter{ 0.9 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{ PbiReadAccuracyFilter{ 0.9, Compare::GREATER_THAN } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,2});
    }
}

TEST(PbiFilterTest, ReadGroupFilterOk)
{
    { // numeric ID
        const auto filter = PbiReadGroupFilter{ -1197849594 };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});

        const auto filter2 = PbiReadGroupFilter{ 200 };
        PbiFilterTests::checkFilterRows(filter2, std::vector<size_t>{});
    }
    { // string ID
        const auto filter = PbiReadGroupFilter{ "b89a4406" };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});

        const auto filter2 = PbiReadGroupFilter{ "b89a4406" };
        PbiFilterTests::checkFilterRows(filter2, std::vector<size_t>{0,1,2,3});
    }
    { // ReadGroupInfo object
        const auto rg = ReadGroupInfo{ "b89a4406" };
        const auto filter = PbiReadGroupFilter{ rg };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    { // multi-ID
        const auto ids = std::vector<int32_t>({-1197849594, 200});
        const auto filter = PbiReadGroupFilter{ ids };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    { // multi-string
        const auto ids = std::vector<std::string>({"b89a4406", "deadbeef"});
        const auto filter = PbiReadGroupFilter{ ids };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    { // multi-ReadGroupInfo
        const auto ids = std::vector<ReadGroupInfo>({ ReadGroupInfo("b89a4406"), ReadGroupInfo("deadbeef")});
        const auto filter = PbiReadGroupFilter{ ids };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
}

TEST(PbiFilterTest, ReferenceEndFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReferenceEndFilter{ 9900 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{ PbiReferenceEndFilter{ 9900, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,3});
    }
}

TEST(PbiFilterTest, ReferenceIdFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ 0 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ 0, Compare::NOT_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto ids = std::vector<int32_t>({0, 42});
        const auto filter = PbiFilter{ PbiReferenceIdFilter{ ids } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
}

TEST(PbiFilterTest, ReferenceNameFilterOk)
{
    const auto bamFile = BamFile{ PbbamTestsConfig::Data_Dir + std::string{ "/group/test2.bam" } };
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
        const auto names = std::vector<std::string>({ "lambda_NEB3011" }); // this file only has 1 :(
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
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{ PbiReferenceStartFilter{ 9200, Compare::GREATER_THAN_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,3});
    }
}

TEST(PbiFilterTest, ZmwFilterOk)
{
    {
        const auto filter = PbiFilter{ PbiZmwFilter{ 14743 } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
    }
    {
        const auto filter = PbiFilter{ PbiZmwFilter{ 14743, Compare::NOT_EQUAL } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto zmws = std::vector<int32_t>({14743,42,200});
        const auto filter = PbiFilter{ PbiZmwFilter{ zmws } };
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0,1,2,3});
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

    for (size_t i = 0; i < PbiFilterTests::shared_index.NumReads(); ++i) {
        EXPECT_EQ(expectedFilter.Accepts(PbiFilterTests::shared_index, i),
                  generatedFilter.Accepts(PbiFilterTests::shared_index, i));
    }
}

TEST(PbiFilterTest, BarcodeListFromDataSetXmlOk)
{
    auto runner = [](const Property& property,
                     const PbiFilter& expectedFilter,
                     const std::vector<size_t>& expectedResults)
    {
        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  expectedResults);
        PbiFilterTests::checkFilterRows(generatedFilter, expectedResults);
    };

    // single barcode
    runner(Property{ "bc", "18", "==" },
           PbiBarcodeFilter{ 18, Compare::EQUAL },
           std::vector<size_t>{1,3});

    // single barcode (bracketed)
    runner(Property{ "bc", "[18]", "==" },
           PbiBarcodeFilter{ 18, Compare::EQUAL },
           std::vector<size_t>{1,3});

    // barcode pair (square brackets)
    runner(Property{ "bc", "[17,18]", "==" },
           PbiBarcodesFilter{ {17, 18}, Compare::EQUAL },
           std::vector<size_t>{1,3});

    // barcode pair (parens)
    runner(Property{ "bc", "(17,18)", "==" },
           PbiBarcodesFilter{ {17, 18}, Compare::EQUAL },
           std::vector<size_t>{1,3});

    // barcode pair (curly brackets)
    runner(Property{ "bc", "{17,18}", "==" },
           PbiBarcodesFilter{ {17, 18}, Compare::EQUAL },
           std::vector<size_t>{1,3});

    // barcode pair (list, but no brackets)
    runner(Property{ "bc", "17,18", "==" },
           PbiBarcodesFilter{ {17, 18}, Compare::EQUAL },
           std::vector<size_t>{1,3});

    // barcode pair - same value
    runner(Property{ "bc", "[18,18]", "==" },
           PbiBarcodesFilter{ {18, 18}, Compare::EQUAL },
           std::vector<size_t>{}); // none share forward & reverse

    auto expectFail = [](const Property& property)
    {
        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        EXPECT_THROW(PbiFilter::FromDataSet(dataset), std::runtime_error);
    };

    // list-ish, but only one value
    expectFail(Property{ "bc", "[18,]", "==" });

    // too many barcodes
    expectFail(Property{ "bc", "[18,18,18]", "==" });
}

TEST(PbiFilterTest, LocalContextFiltersFromDataSetXmlOk)
{
    {   // no adapters or barcodes

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::NO_LOCAL_CONTEXT, Compare::EQUAL };

        // XML: <Property Name="cx" Value="0" Operator="==" />
        Property property("cx", "0", "==");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{0});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{0});
    }
    {   // any adapters or barcodes

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::NO_LOCAL_CONTEXT, Compare::NOT_EQUAL };

        // XML: <Property Name="cx" Value="0" Operator="!=" />
        Property property("cx", "0", "!=");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,2,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,2,3});
    }
    {   // contains adapter_before

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS };

        // XML: <Property Name="cx" Value="1" Operator="&" />
        Property property("cx", "1", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,3});
    }
    {   // contains adapter_before

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS };

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,3});
    }
    {   // contains adapter_after

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER, Compare::CONTAINS };

        // XML: <Property Name="cx" Value="2" Operator="&" />
        Property property("cx", "2", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{2,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{2,3});
    }
    {   // contains adapter_before or adapter_after

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE | LocalContextFlags::ADAPTER_AFTER,
                                       Compare::CONTAINS };

        // XML: <Property Name="cx" Value="3" Operator="&" />
        Property property("cx", "3", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,2,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,2,3});
    }
    {   // contains adapter_before or adapter_after

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE | LocalContextFlags::ADAPTER_AFTER,
                                       Compare::CONTAINS };

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE | ADAPTER_AFTER" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE | ADAPTER_AFTER", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,2,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,2,3});
    }
    {   // contains adapter_before or adapter_after - no whitespace separation

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE | LocalContextFlags::ADAPTER_AFTER,
                                       Compare::CONTAINS };

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE|ADAPTER_AFTER" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE|ADAPTER_AFTER", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,2,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,2,3});
    }
    {   // contains adapter_before or adapter_after - a lot of whitespace separation

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE | LocalContextFlags::ADAPTER_AFTER,
                                       Compare::CONTAINS };

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE        |           ADAPTER_AFTER" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE        |           ADAPTER_AFTER", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,2,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,2,3});
    }
    {   // contains adapter_before or adapter_after, but not both

        const auto expectedFilter = PbiFilter::Union(
        {
            PbiFilter::Intersection(
            {
                PbiLocalContextFilter{ LocalContextFlags::NO_LOCAL_CONTEXT, Compare::NOT_EQUAL },
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS }
            }),
            PbiFilter::Intersection(
            {
                PbiLocalContextFilter{ LocalContextFlags::NO_LOCAL_CONTEXT, Compare::NOT_EQUAL },
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER, Compare::NOT_CONTAINS }
            })
        });

        // XML:
        // <Filters>
        //   <Filter>
        //     <Properties>
        //       <Property Name="cx" Value="0" Operator="!=" />
        //       <Property Name="cx" Value="1" Operator="~" />
        //     </Properties>
        //   </Filter>
        //   <Filter>
        //     <Properties>
        //       <Property Name="cx" Value="0" Operator="!=" />
        //       <Property Name="cx" Value="2" Operator="~" />
        //     </Properties>
        //   </Filter>
        // </Filters>

        auto filter1 = Filter{ };
        filter1.Properties().Add(Property("cx", "0", "!="));
        filter1.Properties().Add(Property("cx", "1", "~"));

        auto filter2 = Filter{ };
        filter2.Properties().Add(Property("cx", "0", "!="));
        filter2.Properties().Add(Property("cx", "2", "~"));

        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter1);
        dataset.Filters().Add(filter2);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,2});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,2});

    }
    {   // contains adapter_before or adapter_after

        const auto expectedFilter = PbiFilter::Union(
        {
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS },
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::CONTAINS }
        });

        // XML:
        // <Filters>
        //   <Filter>
        //     <Properties>
        //       <Property Name="cx" Value="1" Operator="&" />
        //     </Properties>
        //   </Filter>
        //   <Filter>
        //     <Properties>
        //       <Property Name="cx" Value="2" Operator="&" />
        //     </Properties>
        //   </Filter>
        // </Filters>

        auto filter1 = Filter{ };
        filter1.Properties().Add(Property("cx", "1", "&"));

        auto filter2 = Filter{ };
        filter2.Properties().Add(Property("cx", "2", "&"));

        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter1);
        dataset.Filters().Add(filter2);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1,2,3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1,2,3});
    }
    { // adapter_before and adapter_after

        const auto expectedFilter = PbiFilter::Intersection(
        {
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS },
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::CONTAINS }
        });

        // XML:
        // <Property Name="cx" Value="1" Operator="&" />
        // <Property Name="cx" Value="2" Operator="&" />
        Property property1("cx", "1", "&");
        Property property2("cx", "2", "&");

        auto filter = Filter{ };
        filter.Properties().Add(property1);
        filter.Properties().Add(property2);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{3});
    }
    {   // adapter_before, but no adapter_after

        const auto expectedFilter = PbiFilter::Intersection(
        {
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS },
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::NOT_CONTAINS }
        });

        // XML:
        // <Property Name="cx" Value="1" Operator="&" />
        // <Property Name="cx" Value="2" Operator="~" />
        Property property1("cx", "1", "&");
        Property property2("cx", "2", "~");

        auto filter = Filter{ };
        filter.Properties().Add(property1);
        filter.Properties().Add(property2);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{1});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1});
    }
    {   // contains no adapter_before

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS };

        // XML: <Property Name="cx" Value="1" Operator="~" />
        Property property("cx", "1", "~");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{0,2});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{0,2});
    }
    {   // contains no adapter_before or adapter_after

        const auto expectedFilter = PbiFilter::Intersection(
        {
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS },
            PbiLocalContextFilter{ LocalContextFlags::ADAPTER_AFTER,  Compare::NOT_CONTAINS }
        });

        // XML:
        // <Property Name="cx" Value="1" Operator="~" />
        // <Property Name="cx" Value="2" Operator="~" />
        Property property1("cx", "1", "~");
        Property property2("cx", "2", "~");

        auto filter = Filter{ };
        filter.Properties().Add(property1);
        filter.Properties().Add(property2);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{0});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{0});
    }
    {   // contains no adapter_before or adapter_after

        const auto expectedFilter =
                PbiLocalContextFilter{ LocalContextFlags::ADAPTER_BEFORE | LocalContextFlags::ADAPTER_AFTER,
                                       Compare::NOT_CONTAINS };

        // XML: <Property Name="cx" Value="3" Operator="~" />
        Property property("cx", "3", "~");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter,  std::vector<size_t>{0});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{0});
    }
    {   // throws on invalid enum name

        Property property("cx", "DOES_NOT_EXIST", "~");

        auto filter = Filter{ };
        filter.Properties().Add(property);
        DataSet dataset = DataSet{ };
        dataset.Filters().Add(filter);

        EXPECT_THROW(PbiFilter::FromDataSet(dataset), std::runtime_error);
    }
}

// clang-format on
