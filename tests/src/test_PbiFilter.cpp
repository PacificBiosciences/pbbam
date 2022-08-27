#include <pbbam/PbiFilter.h>

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace PbiFilterTests {

// helper structs & methods

static PbiRawData test2Bam_RawIndex()
{
    PbiRawData index;
    index.FileSections(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::REFERENCE);
    index.NumReads(4);

    PbiRawBasicData& subreadData = index.BasicData();
    subreadData.rgId_ = {-1197849594, -1197849594, -1197849594, -1197849594};
    subreadData.qStart_ = {2114, 2579, 4101, 5615};
    subreadData.qEnd_ = {2531, 4055, 5571, 6237};
    subreadData.holeNumber_ = {14743, 14743, 14743, 14743};
    subreadData.readQual_ = {0.901, 0.601, 0.901, 0.601};
    subreadData.ctxtFlag_ = {0, 1, 2, 3};
    subreadData.fileOffset_ = {35651584, 35655125, 35667128, 35679170};

    PbiRawMappedData& mappedData = index.MappedData();
    mappedData.tId_ = {0, 0, 0, 0};
    mappedData.tStart_ = {9507, 8453, 8455, 9291};
    mappedData.tEnd_ = {9903, 9902, 9893, 9900};
    mappedData.aStart_ = {2130, 2581, 4102, 5619};
    mappedData.aEnd_ = {2531, 4055, 5560, 6237};
    mappedData.revStrand_ = {0, 1, 0, 1};
    mappedData.mapQV_ = {254, 254, 254, 254};
    mappedData.nM_ = {384, 1411, 1393, 598};
    mappedData.nMM_ = {0, 0, 0, 0};

    PbiRawReferenceData& referenceData = index.ReferenceData();
    referenceData.entries_.emplace_back(0, 0, 3);
    referenceData.entries_.emplace_back(1);
    referenceData.entries_.emplace_back(PbiReferenceEntry::UNMAPPED_ID);

    return index;
}

static PbiRawData test2Bam_RawBarcodedIndex()
{
    PbiRawData index;
    index.NumReads(4);

    PbiRawBasicData& subreadData = index.BasicData();
    subreadData.rgId_ = {-1197849594, -1197849594, -1197849594, -1197849594};
    subreadData.qStart_ = {2114, 2579, 4101, 5615};
    subreadData.qEnd_ = {2531, 4055, 5571, 6237};
    subreadData.holeNumber_ = {14743, 14743, 14743, 14743};
    subreadData.readQual_ = {0.901, 0.601, 0.901, 0.601};
    subreadData.ctxtFlag_ = {0, 1, 2, 3};
    subreadData.fileOffset_ = {35651584, 35655125, 35667128, 35679170};

    PbiRawMappedData& mappedData = index.MappedData();
    mappedData.tId_ = {0, 0, 0, 0};
    mappedData.tStart_ = {9507, 8453, 8455, 9291};
    mappedData.tEnd_ = {9903, 9902, 9893, 9900};
    mappedData.aStart_ = {2130, 2581, 4102, 5619};
    mappedData.aEnd_ = {2531, 4055, 5560, 6237};
    mappedData.revStrand_ = {0, 1, 0, 1};
    mappedData.mapQV_ = {254, 254, 254, 254};
    mappedData.nM_ = {384, 1411, 1393, 598};
    mappedData.nMM_ = {0, 0, 0, 0};

    PbiRawBarcodeData& barcodeData = index.BarcodeData();
    barcodeData.bcForward_ = {0, 17, 256, 17};
    barcodeData.bcReverse_ = {1, 18, 257, 18};
    barcodeData.bcQual_ = {42, 80, 42, 110};

    PbiRawReferenceData& referenceData = index.ReferenceData();
    referenceData.entries_.emplace_back(0, 0, 3);
    referenceData.entries_.emplace_back(1);
    referenceData.entries_.emplace_back(PbiReferenceEntry::UNMAPPED_ID);

    return index;
}

static const PbiRawData shared_index = test2Bam_RawIndex();
static const PbiRawData shared_barcoded_index = test2Bam_RawBarcodedIndex();

static void checkFilterRows(const PbiFilter& filter, const std::vector<size_t> expectedRows)
{
    if (expectedRows.empty()) {
        for (size_t row = 0; row < shared_index.NumReads(); ++row) {
            EXPECT_FALSE(filter.Accepts(shared_index, row));
        }
    } else {
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(shared_index, row));
        }
    }
}

static void checkFilterBarcodedRows(const PbiFilter& filter, const std::vector<size_t> expectedRows)
{
    if (expectedRows.empty()) {
        for (size_t row = 0; row < shared_barcoded_index.NumReads(); ++row) {
            EXPECT_FALSE(filter.Accepts(shared_barcoded_index, row));
        }
    } else {
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(shared_barcoded_index, row));
        }
    }
}

static void checkFilterInternals(const PbiFilter& filter,
                                 const PbiFilter::CompositionType expectedType,
                                 const size_t expectedNumChildren,
                                 const std::vector<size_t> expectedRows)
{
    EXPECT_EQ(expectedType, filter.Type());
    EXPECT_EQ(expectedNumChildren, filter.NumChildren());
    checkFilterRows(filter, expectedRows);
}

struct SimpleFilter
{
    bool Accepts(const PbiRawData& /* idx */, const size_t /* row */) const
    { /*()idx; ()row;*/
        return true;
    }
};

struct NoncompliantFilter
{};

struct SortUniqueTestFilter
{
    bool Accepts(const PbiRawData& /* idx */, const size_t row) const
    {
        switch (row) {
            case 0:  // fall through
            case 1:  // .
            case 2:  // .
            case 3:  // .
            case 4:  // .
            case 7:  // .
            case 8:
                return true;
            default:
                return false;
        }
    }
};

struct SortUniqueTestFilter2
{
    bool Accepts(const PbiRawData& /* idx */, const size_t row) const
    {
        switch (row) {
            case 3:  // fall through
            case 7:  // .
            case 5:
                return true;
            default:
                return false;
        }
    }
};

static PbiFilter emptyFilter() { return PbiFilter{}; }

static PbiFilter simpleFilter() { return PbiFilter{SimpleFilter{}}; }

}  // namespace PbiFilterTests

TEST(BAM_PbiFilter, default_filter_accepts_all)
{
    const PbiFilter filter;
    PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 0,
                                         std::vector<size_t>{0, 1, 2, 3});
}

TEST(BAM_PbiFilter, can_compose_with_child_filters)
{
    PbiFilter filter;
    filter.Add(PbiFilter{});
    PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1,
                                         std::vector<size_t>{0, 1, 2, 3});
}

TEST(BAM_PbiFilter, can_add_user_defined_filter_types)
{
    {  // ctor
        PbiFilter filter{PbiFilterTests::SimpleFilter{}};
        PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
    }
    {  // Add
        PbiFilter filter;
        filter.Add(PbiFilterTests::SimpleFilter{});
        PbiFilterTests::checkFilterInternals(filter, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
    }

    //    PbiFilter shouldNotCompile = PbiFilter{ PbiFilterTests::NoncompliantFilter{ } };                       // <-- when uncommented, should not compile
    //    PbiFilter shouldNotCompileEither; shouldNotCompileEither.Add(PbiFilterTests::NoncompliantFilter{ });   // <-- when uncommented, should not compile
}

TEST(BAM_PbiFilter, copied_filters_are_equivalent)
{
    {  // empty
        const PbiFilter original;

        const PbiFilter copyCtor{original};
        PbiFilter copyAssign;
        copyAssign = original;

        PbiFilterTests::checkFilterInternals(original, PbiFilter::INTERSECT, 0,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(copyCtor, PbiFilter::INTERSECT, 0,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 0,
                                             std::vector<size_t>{0, 1, 2, 3});
    }
    {  // with children
        const PbiFilter original{PbiFilterTests::SimpleFilter{}};

        PbiFilter copyCtor{original};
        PbiFilter copyAssign;
        copyAssign = original;

        PbiFilterTests::checkFilterInternals(original, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(copyCtor, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(copyAssign, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
    }
}

TEST(BAM_PbiFilter, moved_filters_are_equivalent)
{
    {  // empty
        const PbiFilter original;

        PbiFilter moveCtor(PbiFilterTests::emptyFilter());
        PbiFilter moveAssign;
        moveAssign = PbiFilterTests::emptyFilter();

        PbiFilterTests::checkFilterInternals(original, PbiFilter::INTERSECT, 0,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(moveCtor, PbiFilter::INTERSECT, 0,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 0,
                                             std::vector<size_t>{0, 1, 2, 3});
    }
    {  // with children
        const auto original = PbiFilterTests::simpleFilter();

        PbiFilter moveCtor(PbiFilterTests::simpleFilter());
        PbiFilter moveAssign;
        moveAssign = PbiFilterTests::simpleFilter();

        PbiFilterTests::checkFilterInternals(original, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(moveCtor, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
        PbiFilterTests::checkFilterInternals(moveAssign, PbiFilter::INTERSECT, 1,
                                             std::vector<size_t>{0, 1, 2, 3});
    }
}

TEST(BAM_PbiFilter, SortsAndUniquesChildFilterResultsOk)
{
    const auto childFilter = PbiFilterTests::SortUniqueTestFilter{};
    const auto filter = PbiFilter{childFilter};
    PbiFilterTests::checkFilterRows(childFilter, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
    PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3, 4, 7, 8});
}

TEST(BAM_PbiFilter, can_filter_union_of_child_filters)
{
    { // empty
     {// copy
      const auto emptyFilter = PbiFilterTests::emptyFilter();
    const auto emptyFilter2 = PbiFilterTests::emptyFilter();
    const auto u = PbiFilter::Union({emptyFilter, emptyFilter2});
    PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0, 1, 2, 3});
}
{  // move
    const auto u = PbiFilter::Union({PbiFilter{}, PbiFilter{}});
    PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0, 1, 2, 3});
}
}

{ // with (no-data) children - just checking composition
 {// copy
  const auto simpleFilter = PbiFilterTests::SimpleFilter{};
const auto simpleFilter2 = PbiFilterTests::SimpleFilter{};
const auto u = PbiFilter::Union({simpleFilter, simpleFilter2});
PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0, 1, 2, 3});
}
{  // move
    const auto u =
        PbiFilter::Union({PbiFilterTests::SimpleFilter{}, PbiFilterTests::SimpleFilter{}});
    PbiFilterTests::checkFilterInternals(u, PbiFilter::UNION, 2, std::vector<size_t>{0, 1, 2, 3});
}
}

{  // 2-child union, results sorted & unique-d by PbiFilter

    const auto child1 = PbiFilterTests::SortUniqueTestFilter{};
    const auto child2 = PbiFilterTests::SortUniqueTestFilter2{};
    const auto u = PbiFilter::Union({child1, child2});

    PbiFilterTests::checkFilterRows(child1, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
    PbiFilterTests::checkFilterRows(child2, std::vector<size_t>{3, 7, 5});
    PbiFilterTests::checkFilterRows(u, std::vector<size_t>{0, 1, 2, 3, 4, 5, 7, 8});
}
}

TEST(BAM_PbiFilter, can_filter_intersection_of_child_filters)
{
    { // empty
     {// copy
      const auto emptyFilter = PbiFilterTests::emptyFilter();
    const auto emptyFilter2 = PbiFilterTests::emptyFilter();
    const auto i = PbiFilter::Intersection({emptyFilter, emptyFilter2});
    PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2,
                                         std::vector<size_t>{0, 1, 2, 3});
}
{  // move
    const auto i = PbiFilter::Intersection({PbiFilter{}, PbiFilter{}});
    PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2,
                                         std::vector<size_t>{0, 1, 2, 3});
}
}

{ // with (no-data) children - just checking composition
 {// copy
  const auto simpleFilter = PbiFilterTests::SimpleFilter{};
const auto simpleFilter2 = PbiFilterTests::SimpleFilter{};
const auto i = PbiFilter::Intersection({simpleFilter, simpleFilter2});
PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2, std::vector<size_t>{0, 1, 2, 3});
}
{  // move
    const auto i =
        PbiFilter::Intersection({PbiFilterTests::SimpleFilter{}, PbiFilterTests::SimpleFilter{}});
    PbiFilterTests::checkFilterInternals(i, PbiFilter::INTERSECT, 2,
                                         std::vector<size_t>{0, 1, 2, 3});
}
}

{  // 2-child intersect, sorted & unique-d by PbiFilter

    const auto child1 = PbiFilterTests::SortUniqueTestFilter{};
    const auto child2 = PbiFilterTests::SortUniqueTestFilter2{};
    const auto i = PbiFilter::Intersection({child1, child2});

    PbiFilterTests::checkFilterRows(child1, std::vector<size_t>{2, 7, 0, 3, 4, 1, 8});
    PbiFilterTests::checkFilterRows(child2, std::vector<size_t>{3, 7, 5});
    PbiFilterTests::checkFilterRows(i, std::vector<size_t>{3, 7});
}
}

TEST(BAM_PbiFilter, can_filter_on_aligned_end)
{
    {
        const auto filter = PbiFilter{PbiAlignedEndFilter{4055}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{PbiAlignedEndFilter{4055, Compare::NOT_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 2, 3});
    }
    {
        const auto filter = PbiFilter{PbiAlignedEndFilter{4000, Compare::LESS_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
    {
        const auto filter = PbiFilter{PbiAlignedEndFilter{5560, Compare::GREATER_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{PbiAlignedEndFilter{5560, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2, 3});
    }

    {
        const auto filter = PbiFilter{PbiAlignedEndFilter{7000, Compare::GREATER_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(BAM_PbiFilter, can_filter_on_aligned_length)
{
    {
        const auto filter = PbiFilter{PbiAlignedLengthFilter{500, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2, 3});
    }
    {
        const auto filter = PbiFilter{PbiAlignedLengthFilter{1000, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2});
    }
}

TEST(BAM_PbiFilter, can_filter_on_aligned_start)
{
    {
        const auto filter = PbiFilter{PbiAlignedStartFilter{2600, Compare::LESS_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1});
    }
    {
        const auto filter = PbiFilter{PbiAlignedStartFilter{4102, Compare::GREATER_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{PbiAlignedStartFilter{4102, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2, 3});
    }
    {
        const auto filter = PbiFilter{PbiAlignedStartFilter{6000, Compare::GREATER_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(BAM_PbiFilter, can_filter_on_aligned_strand)
{
    {
        const auto filter = PbiFilter{PbiAlignedStrandFilter{Data::Strand::FORWARD}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 2});
    }
    {
        const auto filter = PbiFilter{PbiAlignedStrandFilter{Data::Strand::REVERSE}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 3});
    }
    {
        const auto filter = PbiFilter{PbiAlignedStrandFilter{
            Data::Strand::FORWARD, Compare::NOT_EQUAL}};  // same as Data::Strand::REVERSE
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 3});
    }

    // unsupported compare types throw
    EXPECT_THROW(PbiAlignedStrandFilter(Data::Strand::FORWARD, Compare::LESS_THAN),
                 std::runtime_error);
    EXPECT_THROW(PbiAlignedStrandFilter(Data::Strand::FORWARD, Compare::LESS_THAN_EQUAL),
                 std::runtime_error);
    EXPECT_THROW(PbiAlignedStrandFilter(Data::Strand::FORWARD, Compare::GREATER_THAN),
                 std::runtime_error);
    EXPECT_THROW(PbiAlignedStrandFilter(Data::Strand::FORWARD, Compare::GREATER_THAN_EQUAL),
                 std::runtime_error);
}

TEST(BAM_PbiFilter, can_filter_on_single_barcode)
{
    {
        const auto filter = PbiFilter{PbiBarcodeFilter{17}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
    {
        const auto filter = PbiFilter{PbiBarcodeFilter{18}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
    {
        const auto filter = PbiFilter{PbiBarcodeFilter{0}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{0});
    }
}

TEST(BAM_PbiFilter, can_filter_on_barcode_forward)
{
    {
        const auto filter = PbiFilter{PbiBarcodeForwardFilter{17}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
    {
        const auto filter = PbiFilter{PbiBarcodeForwardFilter{400}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{PbiBarcodeForwardFilter{{0, 256}}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{0, 2});
    }
    {
        //blacklist
        const auto filter = PbiFilter{PbiBarcodeForwardFilter{{0, 256}, Compare::NOT_CONTAINS}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_barcode_quality)
{
    {
        const auto filter = PbiFilter{PbiBarcodeQualityFilter{80, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
    {
        const auto filter = PbiFilter{PbiBarcodeQualityFilter{40, Compare::LESS_THAN}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{});
    }
}

TEST(BAM_PbiFilter, can_filter_on_barcode_reverse)
{
    {
        const auto filter = PbiFilter{PbiBarcodeReverseFilter{18}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
    {
        const auto filter = PbiFilter{PbiBarcodeReverseFilter{400}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{PbiBarcodeReverseFilter{{1, 257}}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{0, 2});
    }
    {
        // blacklist
        const auto filter = PbiFilter{PbiBarcodeReverseFilter{{1, 257}, Compare::NOT_CONTAINS}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_barcode_pair)
{
    {
        const auto filter = PbiFilter{PbiBarcodesFilter{17, 18}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
    {
        const auto filter = PbiFilter{PbiBarcodesFilter{17, 19}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{PbiBarcodesFilter{std::make_pair(17, 18)}};
        PbiFilterTests::checkFilterBarcodedRows(filter, std::vector<size_t>{1, 3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_identity)
{
    {
        const auto filter = PbiFilter{PbiIdentityFilter{0.95, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_local_context_flags)
{
    {  // == NO_LOCAL_CONTEXT
        const auto filter =
            PbiFilter{PbiLocalContextFilter{Data::LocalContextFlags::NO_LOCAL_CONTEXT}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
    {  // != ADAPTER_BEFORE (exact match)
        const auto filter = PbiFilter{
            PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 2, 3});
    }
    {  // contains ADAPTER_BEFORE
        const auto filter = PbiFilter{
            PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 3});
    }
    {  // does not contain ADAPTER_BEFORE
        const auto filter = PbiFilter{
            PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 2});
    }
    {  // include both ADAPTER_BEFORE and ADAPTER_AFTER
        const auto filter = PbiFilter::Intersection(
            {PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS},
             PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::CONTAINS}});
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {  // exclude both ADAPTER_BEFORE and ADAPTER_AFTER
        const auto filter = PbiFilter::Intersection(
            {PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS},
             PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::NOT_CONTAINS}});
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
    {  // include everything with either ADAPTER_BEFORE or ADAPTER_AFTER
        const auto filter = PbiFilter::Union(
            {PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS},
             PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::CONTAINS}});
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2, 3});
    }
    {  // include everything with either ADAPTER_BEFORE or ADAPTER_AFTER, but not both
        const auto filter = PbiFilter::Intersection(
            {PbiLocalContextFilter{Data::LocalContextFlags::NO_LOCAL_CONTEXT, Compare::NOT_EQUAL},
             PbiFilter::Union({PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE,
                                                     Compare::NOT_CONTAINS},
                               PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER,
                                                     Compare::NOT_CONTAINS}})});
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2});
    }
}

TEST(BAM_PbiFilter, can_filter_on_map_quality)
{
    {
        const auto filter = PbiFilter{PbiMapQualityFilter{254}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {
        const auto filter = PbiFilter{PbiMapQualityFilter{254, Compare::NOT_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(BAM_PbiFilter, can_filter_on_movie_name)
{
    const auto bamFile = BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}};
    const auto index = PbiRawData{bamFile.PacBioIndexFilename()};

    {
        const auto filter = PbiFilter{
            PbiMovieNameFilter{"m140905_042212_sidney_c100564852550000001823085912221377_s1_X0"}};
        const auto expectedRows = std::vector<size_t>{0, 1, 2, 3};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
    {
        const auto filter = PbiFilter{PbiMovieNameFilter{"does_not_exist"}};
        const auto expectedRows = std::vector<size_t>{};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
    {
        const auto names = std::vector<std::string>{
            "does_not_exist", "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0"};
        const auto filter = PbiFilter{PbiMovieNameFilter{names}};
        const auto expectedRows = std::vector<size_t>{0, 1, 2, 3};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
    {
        // blacklist
        const auto names = std::vector<std::string>{
            "does_not_exist", "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0"};
        const auto filter = PbiFilter{PbiMovieNameFilter{names, Compare::NOT_CONTAINS}};
        for (size_t row = 0; row < index.NumReads(); ++row) {
            EXPECT_FALSE(filter.Accepts(index, row));
        }
    }
}

TEST(BAM_PbiFilter, can_filter_on_num_deleted_bases)
{
    // del: { 12, 38, 45, 11} - calculated from raw data, not stored directly in testing object or read from PBI file

    {
        const auto filter = PbiFilter{PbiNumDeletedBasesFilter{12, Compare::LESS_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 3});
    }
    {
        const auto filter = PbiFilter{PbiNumDeletedBasesFilter{45, Compare::EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2});
    }
}

TEST(BAM_PbiFilter, can_filter_on_num_inserted_bases)
{
    // ins: { 17, 63, 65, 20 }  - calculated from raw data, not stored directly testing object or read from PBI file

    {
        const auto filter = PbiFilter{PbiNumInsertedBasesFilter{63, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2});
    }
    {
        const auto filter = PbiFilter{PbiNumInsertedBasesFilter{17, Compare::NOT_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2, 3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_num_matches)
{
    {
        const auto filter = PbiFilter{PbiNumMatchesFilter{1000, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2});
    }
    {
        const auto filter = PbiFilter{PbiNumMatchesFilter{400, Compare::LESS_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0});
    }
}

TEST(BAM_PbiFilter, can_filter_on_num_mismatches)
{
    {
        const auto filter = PbiFilter{PbiNumMismatchesFilter{0, Compare::EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {
        const auto filter = PbiFilter{PbiNumMismatchesFilter{0, Compare::NOT_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(BAM_PbiFilter, can_filter_on_num_subreads)
{
    PbiRawData index;
    index.NumReads(21);

    PbiRawBasicData& subreadData = index.BasicData();
    subreadData.rgId_ = std::vector<int32_t>(21, 0);
    subreadData.qStart_ = std::vector<int32_t>(21, 0);
    subreadData.qEnd_ = std::vector<int32_t>(21, 0);
    subreadData.readQual_ = std::vector<float>(21, 0);
    subreadData.ctxtFlag_ = std::vector<uint8_t>(21, 0);
    subreadData.fileOffset_ = std::vector<int64_t>(21, 0);
    subreadData.fileNumber_ = std::vector<uint16_t>(21, 0);
    subreadData.holeNumber_ = {0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 6};

    {  // >= 3
        const auto filter = PbiNumSubreadsFilter{3, Compare::GREATER_THAN_EQUAL};
        const std::vector<size_t> expectedRows{0, 1, 2, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
    {  // < 3
        const auto filter = PbiNumSubreadsFilter{3, Compare::LESS_THAN};
        const std::vector<size_t> expectedRows{3, 4, 10, 18, 19, 20};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
    {  // == 1
        const auto filter = PbiNumSubreadsFilter{1, Compare::EQUAL};
        const std::vector<size_t> expectedRows{10, 20};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
}

TEST(BAM_PbiFilter, can_filter_on_query_end)
{
    {
        const auto filter = PbiFilter{PbiQueryEndFilter{4055}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{PbiQueryEndFilter{6200, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_query_length)
{
    {
        const auto filter = PbiFilter{PbiQueryLengthFilter{500, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2, 3});
    }
    {
        const auto filter = PbiFilter{PbiQueryLengthFilter{1000, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 2});
    }
}

TEST(BAM_PbiFilter, can_filter_on_query_name)
{
    {
        const auto filter = PbiFilter{PbiQueryNameFilter{
            "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055"}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{PbiQueryNameFilter{
            "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }

    {
        const auto filter = PbiFilter{PbiQueryNameFilter{"does_not_exist/0/0_0"}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto names = std::vector<std::string>{
            "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055",
            "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"};
        const auto filter = PbiFilter{PbiQueryNameFilter{names}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1, 3});
    }
}

TEST(BAM_PbiFilter, throws_on_invalid_on_query_name)
{
    // invalid QNAME syntax throws
    EXPECT_THROW(
        {
            const auto filter = PbiFilter{PbiQueryNameFilter{""}};
            PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            const auto filter = PbiFilter{PbiQueryNameFilter{"foo"}};
            PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            const auto filter = PbiFilter{PbiQueryNameFilter{"foo/bar"}};
            PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            const auto filter = PbiFilter{PbiQueryNameFilter{"foo/bar/baz_bam"}};
            PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
        },
        std::exception);  // come back to see why this is not runtime_error but something else
}

TEST(BAM_PbiFilter, can_filter_on_query_start)
{
    {
        const auto filter = PbiFilter{PbiQueryStartFilter{4101}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{2});
    }
    {
        const auto filter = PbiFilter{PbiQueryStartFilter{5000}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{PbiQueryStartFilter{5000, Compare::GREATER_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_read_accuracy)
{
    {
        const auto filter = PbiFilter{PbiReadAccuracyFilter{0.9}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto filter = PbiFilter{PbiReadAccuracyFilter{0.9, Compare::GREATER_THAN}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 2});
    }
}

TEST(BAM_PbiFilter, can_filter_on_read_group)
{
    {  // numeric ID
        const auto filter = PbiReadGroupFilter{-1197849594};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});

        const auto filter2 = PbiReadGroupFilter{200};
        PbiFilterTests::checkFilterRows(filter2, std::vector<size_t>{});
    }
    {  // string ID
        const auto filter = PbiReadGroupFilter{"b89a4406"};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});

        const auto filter2 = PbiReadGroupFilter{"b89a4406"};
        PbiFilterTests::checkFilterRows(filter2, std::vector<size_t>{0, 1, 2, 3});
    }
    {  // ReadGroupInfo object
        const auto rg = ReadGroupInfo{"b89a4406"};
        const auto filter = PbiReadGroupFilter{rg};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {  // multi-ID
        const auto ids = std::vector<int32_t>({-1197849594, 200});
        const auto filter = PbiReadGroupFilter{ids};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {  // multi-ID blacklist
        const auto ids = std::vector<int32_t>({-1197849594, 200});
        const auto filter = PbiReadGroupFilter{ids, Compare::NOT_CONTAINS};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {  // multi-string
        const auto ids = std::vector<std::string>({"b89a4406", "deadbeef"});
        const auto filter = PbiReadGroupFilter{ids};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {  // multi-ReadGroupInfo
        const auto ids =
            std::vector<ReadGroupInfo>({ReadGroupInfo("b89a4406"), ReadGroupInfo("deadbeef")});
        const auto filter = PbiReadGroupFilter{ids};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_reference_end)
{
    {
        const auto filter = PbiFilter{PbiReferenceEndFilter{9900}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{3});
    }
    {
        const auto filter = PbiFilter{PbiReferenceEndFilter{9900, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_reference_id)
{
    {
        const auto filter = PbiFilter{PbiReferenceIdFilter{0}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {
        const auto filter = PbiFilter{PbiReferenceIdFilter{0, Compare::NOT_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto ids = std::vector<int32_t>({0, 42});
        const auto filter = PbiFilter{PbiReferenceIdFilter{ids}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {
        const auto ids = std::vector<int32_t>({0});
        const auto filter = PbiFilter{PbiReferenceIdFilter{ids, Compare::NOT_CONTAINS}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(BAM_PbiFilter, can_filter_on_reference_name)
{
    const auto bamFile = BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}};
    const auto index = PbiRawData{bamFile.PacBioIndexFilename()};

    {
        const auto filter = PbiFilter{PbiReferenceNameFilter{"lambda_NEB3011"}};
        const auto expectedRows = std::vector<size_t>{0, 1, 2, 3};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
    {
        const auto filter = PbiFilter{PbiReferenceNameFilter{"lambda_NEB3011", Compare::NOT_EQUAL}};
        const auto expectedRows = std::vector<size_t>{};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }
    {
        const auto names = std::vector<std::string>({"lambda_NEB3011"});  // this file only has 1 :(
        const auto filter = PbiFilter{PbiReferenceNameFilter{names}};
        const auto expectedRows = std::vector<size_t>{0, 1, 2, 3};
        for (size_t row : expectedRows) {
            EXPECT_TRUE(filter.Accepts(index, row));
        }
    }

    // unsupported compare types throw
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::LESS_THAN), std::runtime_error);
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::LESS_THAN_EQUAL), std::runtime_error);
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::GREATER_THAN), std::runtime_error);
    EXPECT_THROW(PbiReferenceNameFilter("foo", Compare::GREATER_THAN_EQUAL), std::runtime_error);
}

TEST(BAM_PbiFilter, can_filter_on_reference_start)
{
    {
        const auto filter = PbiFilter{PbiReferenceStartFilter{8453}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{1});
    }
    {
        const auto filter = PbiFilter{PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 3});
    }
}

TEST(BAM_PbiFilter, can_filter_on_zmw)
{
    {
        const auto filter = PbiFilter{PbiZmwFilter{14743}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {
        // blacklist
        const auto filter = PbiFilter{PbiZmwFilter{14743, Compare::NOT_EQUAL}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
    {
        const auto zmws = std::vector<int32_t>({14743, 42, 200});
        const auto filter = PbiFilter{PbiZmwFilter{zmws}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{0, 1, 2, 3});
    }
    {
        // blacklist
        const auto zmws = std::vector<int32_t>{14743};
        const auto filter = PbiFilter{PbiZmwFilter{zmws, Compare::NOT_CONTAINS}};
        PbiFilterTests::checkFilterRows(filter, std::vector<size_t>{});
    }
}

TEST(BAM_PbiFilter, can_load_from_dataset)
{
    const auto expectedFilter = PbiFilter::Union(
        {PbiFilter::Intersection(
             {PbiZmwFilter{14743}, PbiReadAccuracyFilter{0.9, Compare::GREATER_THAN_EQUAL}}),

         PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}});

    auto properties1 = Properties{};
    properties1.Add(Property{"zm", "14743", "=="});
    properties1.Add(Property{"rq", "0.9", ">="});

    auto datasetFilter1 = Filter{};
    datasetFilter1.Properties(properties1);

    auto properties2 = Properties{};
    properties2.Add(Property{"pos", "9200", ">="});

    auto datasetFilter2 = Filter{};
    datasetFilter2.Properties(properties2);

    auto datasetFilters = Filters{};
    datasetFilters.Add(datasetFilter1);
    datasetFilters.Add(datasetFilter2);
    auto dataset = DataSet{};
    dataset.Filters(datasetFilters);

    const auto generatedFilter = PbiFilter::FromDataSet(dataset);

    for (size_t i = 0; i < PbiFilterTests::shared_index.NumReads(); ++i) {
        EXPECT_EQ(expectedFilter.Accepts(PbiFilterTests::shared_index, i),
                  generatedFilter.Accepts(PbiFilterTests::shared_index, i));
    }
}

TEST(BAM_PbiFilter, can_load_from_dataset_with_barcode_list)
{
    auto runner = [](const Property& property, const PbiFilter& expectedFilter,
                     const std::vector<size_t>& expectedResults) {
        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterBarcodedRows(expectedFilter, expectedResults);
        PbiFilterTests::checkFilterBarcodedRows(generatedFilter, expectedResults);
    };

    // single barcode
    runner(Property{"bc", "18", "=="}, PbiBarcodeFilter{18, Compare::EQUAL},
           std::vector<size_t>{1, 3});

    // single barcode (bracketed)
    runner(Property{"bc", "[18]", "=="}, PbiBarcodeFilter{18, Compare::EQUAL},
           std::vector<size_t>{1, 3});

    // barcode pair (square brackets)
    runner(Property{"bc", "[17,18]", "=="}, PbiBarcodesFilter{{17, 18}, Compare::EQUAL},
           std::vector<size_t>{1, 3});

    // barcode pair (parens)
    runner(Property{"bc", "(17,18)", "=="}, PbiBarcodesFilter{{17, 18}, Compare::EQUAL},
           std::vector<size_t>{1, 3});

    // barcode pair (curly brackets)
    runner(Property{"bc", "{17,18}", "=="}, PbiBarcodesFilter{{17, 18}, Compare::EQUAL},
           std::vector<size_t>{1, 3});

    // barcode pair (list, but no brackets)
    runner(Property{"bc", "17,18", "=="}, PbiBarcodesFilter{{17, 18}, Compare::EQUAL},
           std::vector<size_t>{1, 3});

    // barcode pair - same value
    runner(Property{"bc", "[18,18]", "=="}, PbiBarcodesFilter{{18, 18}, Compare::EQUAL},
           std::vector<size_t>{});  // none share forward & reverse

    auto expectFail = [](const Property& property) {
        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        EXPECT_THROW(PbiFilter::FromDataSet(dataset), std::runtime_error);
    };

    // list-ish, but only one value
    expectFail(Property{"bc", "[18,]", "=="});

    // too many barcodes
    expectFail(Property{"bc", "[18,18,18]", "=="});
}

TEST(BAM_PbiFilter, can_load_from_dataset_with_local_context)
{
    {  // no adapters or barcodes

        const auto expectedFilter =
            PbiLocalContextFilter{Data::LocalContextFlags::NO_LOCAL_CONTEXT, Compare::EQUAL};

        // XML: <Property Name="cx" Value="0" Operator="==" />
        Property property("cx", "0", "==");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterBarcodedRows(expectedFilter, std::vector<size_t>{0});
        PbiFilterTests::checkFilterBarcodedRows(generatedFilter, std::vector<size_t>{0});
    }
    {  // any adapters or barcodes

        const auto expectedFilter =
            PbiLocalContextFilter{Data::LocalContextFlags::NO_LOCAL_CONTEXT, Compare::NOT_EQUAL};

        // XML: <Property Name="cx" Value="0" Operator="!=" />
        Property property("cx", "0", "!=");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 2, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 2, 3});
    }
    {  // contains adapter_before

        const auto expectedFilter =
            PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS};

        // XML: <Property Name="cx" Value="1" Operator="&" />
        Property property("cx", "1", "&");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 3});
    }
    {  // contains adapter_before

        const auto expectedFilter =
            PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS};

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE", "&");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 3});
    }
    {  // contains adapter_after

        const auto expectedFilter =
            PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::CONTAINS};

        // XML: <Property Name="cx" Value="2" Operator="&" />
        Property property("cx", "2", "&");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{2, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{2, 3});
    }
    {  // contains adapter_before or adapter_after

        const auto expectedFilter = PbiLocalContextFilter{
            Data::LocalContextFlags::ADAPTER_BEFORE | Data::LocalContextFlags::ADAPTER_AFTER,
            Compare::CONTAINS};

        // XML: <Property Name="cx" Value="3" Operator="&" />
        Property property("cx", "3", "&");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 2, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 2, 3});
    }
    {  // contains adapter_before or adapter_after

        const auto expectedFilter = PbiLocalContextFilter{
            Data::LocalContextFlags::ADAPTER_BEFORE | Data::LocalContextFlags::ADAPTER_AFTER,
            Compare::CONTAINS};

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE | ADAPTER_AFTER" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE | ADAPTER_AFTER", "&");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 2, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 2, 3});
    }
    {  // contains adapter_before or adapter_after - no whitespace separation

        const auto expectedFilter = PbiLocalContextFilter{
            Data::LocalContextFlags::ADAPTER_BEFORE | Data::LocalContextFlags::ADAPTER_AFTER,
            Compare::CONTAINS};

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE|ADAPTER_AFTER" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE|ADAPTER_AFTER", "&");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 2, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 2, 3});
    }
    {  // contains adapter_before or adapter_after - a lot of whitespace separation

        const auto expectedFilter = PbiLocalContextFilter{
            Data::LocalContextFlags::ADAPTER_BEFORE | Data::LocalContextFlags::ADAPTER_AFTER,
            Compare::CONTAINS};

        // XML: <Property Name="cx" Value="ADAPTER_BEFORE        |           ADAPTER_AFTER" Operator="&" />
        Property property("cx", "ADAPTER_BEFORE        |           ADAPTER_AFTER", "&");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 2, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 2, 3});
    }
    {  // contains adapter_before or adapter_after, but not both

        const auto expectedFilter =
            PbiFilter::Union({PbiFilter::Intersection(
                                  {PbiLocalContextFilter{Data::LocalContextFlags::NO_LOCAL_CONTEXT,
                                                         Compare::NOT_EQUAL},
                                   PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE,
                                                         Compare::NOT_CONTAINS}}),
                              PbiFilter::Intersection(
                                  {PbiLocalContextFilter{Data::LocalContextFlags::NO_LOCAL_CONTEXT,
                                                         Compare::NOT_EQUAL},
                                   PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER,
                                                         Compare::NOT_CONTAINS}})});

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

        auto filter1 = Filter{};
        filter1.Properties().Add(Property("cx", "0", "!="));
        filter1.Properties().Add(Property("cx", "1", "~"));

        auto filter2 = Filter{};
        filter2.Properties().Add(Property("cx", "0", "!="));
        filter2.Properties().Add(Property("cx", "2", "~"));

        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter1);
        dataset.Filters().Add(filter2);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 2});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 2});
    }
    {  // contains adapter_before or adapter_after

        const auto expectedFilter = PbiFilter::Union(
            {PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS},
             PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::CONTAINS}});

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

        auto filter1 = Filter{};
        filter1.Properties().Add(Property("cx", "1", "&"));

        auto filter2 = Filter{};
        filter2.Properties().Add(Property("cx", "2", "&"));

        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter1);
        dataset.Filters().Add(filter2);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1, 2, 3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1, 2, 3});
    }
    {  // adapter_before and adapter_after

        const auto expectedFilter = PbiFilter::Intersection(
            {PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS},
             PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::CONTAINS}});

        // XML:
        // <Property Name="cx" Value="1" Operator="&" />
        // <Property Name="cx" Value="2" Operator="&" />
        Property property1("cx", "1", "&");
        Property property2("cx", "2", "&");

        auto filter = Filter{};
        filter.Properties().Add(property1);
        filter.Properties().Add(property2);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{3});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{3});
    }
    {  // adapter_before, but no adapter_after

        const auto expectedFilter = PbiFilter::Intersection(
            {PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::CONTAINS},
             PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::NOT_CONTAINS}});

        // XML:
        // <Property Name="cx" Value="1" Operator="&" />
        // <Property Name="cx" Value="2" Operator="~" />
        Property property1("cx", "1", "&");
        Property property2("cx", "2", "~");

        auto filter = Filter{};
        filter.Properties().Add(property1);
        filter.Properties().Add(property2);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{1});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{1});
    }
    {  // contains no adapter_before

        const auto expectedFilter =
            PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS};

        // XML: <Property Name="cx" Value="1" Operator="~" />
        Property property("cx", "1", "~");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{0, 2});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{0, 2});
    }
    {  // contains no adapter_before or adapter_after

        const auto expectedFilter = PbiFilter::Intersection(
            {PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_BEFORE, Compare::NOT_CONTAINS},
             PbiLocalContextFilter{Data::LocalContextFlags::ADAPTER_AFTER, Compare::NOT_CONTAINS}});

        // XML:
        // <Property Name="cx" Value="1" Operator="~" />
        // <Property Name="cx" Value="2" Operator="~" />
        Property property1("cx", "1", "~");
        Property property2("cx", "2", "~");

        auto filter = Filter{};
        filter.Properties().Add(property1);
        filter.Properties().Add(property2);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{0});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{0});
    }
    {  // contains no adapter_before or adapter_after

        const auto expectedFilter = PbiLocalContextFilter{
            Data::LocalContextFlags::ADAPTER_BEFORE | Data::LocalContextFlags::ADAPTER_AFTER,
            Compare::NOT_CONTAINS};

        // XML: <Property Name="cx" Value="3" Operator="~" />
        Property property("cx", "3", "~");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        const auto generatedFilter = PbiFilter::FromDataSet(dataset);
        PbiFilterTests::checkFilterRows(expectedFilter, std::vector<size_t>{0});
        PbiFilterTests::checkFilterRows(generatedFilter, std::vector<size_t>{0});
    }
    {  // throws on invalid enum name

        Property property("cx", "DOES_NOT_EXIST", "~");

        auto filter = Filter{};
        filter.Properties().Add(property);
        DataSet dataset = DataSet{};
        dataset.Filters().Add(filter);

        EXPECT_THROW(PbiFilter::FromDataSet(dataset), std::runtime_error);
    }
}
