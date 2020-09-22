// Author: David Seifert

#include <pbbam/DataSet.h>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace GenomicIntervalsTests {

const std::string inputDir{PbbamTestsConfig::Data_Dir + "/test_GenomicIntervals/"};

}  // namespace GenomicIntervalsTests

TEST(BAM_DataSetGenomicIntervals, fetches_intervals_with_no_filter)
{
    // vanilla AlignmentSet, no filters
    DataSet ds{GenomicIntervalsTests::inputDir + "no_filter.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 0, 20}, {"contig2", 0, 10}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, fetches_no_intervals_with_empty_input)
{
    // interval contig1:[5, 5), i.e., empty, yet both offsets are within range
    DataSet ds{GenomicIntervalsTests::inputDir + "empty.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct;
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, fetches_no_intervals_if_out_of_range)
{
    // interval contig1:[1000, 10000), i.e., empty, as the selected range
    // lies above the contig1 size of 20
    DataSet ds{GenomicIntervalsTests::inputDir + "out_of_range.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct;
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_single_normal_interval)
{
    // interval contig1:[3, 10)
    DataSet ds{GenomicIntervalsTests::inputDir + "single_interval.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 10}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_whole_contig_with_integers)
{
    // interval contig1:[0, 20), i.e., select the whole contig
    DataSet ds{GenomicIntervalsTests::inputDir + "whole_contig.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 0, 20}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_whole_contig_with_name_only)
{
    // interval contig1, i.e., select the whole contig, without a range filter
    DataSet ds{GenomicIntervalsTests::inputDir + "contig_name_only.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 0, 20}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_single_interval_less_than_equal)
{
    // interval contig1:[3, 11), test "tstart <=" relation
    DataSet ds{GenomicIntervalsTests::inputDir + "single_interval_start_lte.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 11}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_single_interval_greater_than_equal)
{
    // interval contig1:[2, 10), test "tend >=" relation
    DataSet ds{GenomicIntervalsTests::inputDir + "single_interval_end_gte.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 2, 10}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_disjoint_intervals)
{
    // interval contig1:[3, 7),[13, 17), test that disjoint intervals remain disjoint
    DataSet ds{GenomicIntervalsTests::inputDir + "disjoint_intervals.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 7}, {"contig1", 13, 17}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_adjacent_intervals)
{
    // interval contig1:[3, 17), test that intervals [3, 10) and [10, 17)
    // get merged into a single overall interval
    DataSet ds{GenomicIntervalsTests::inputDir + "adjacent_intervals.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 17}};
    EXPECT_EQ(correct, result);
}

TEST(BAM_DataSetGenomicIntervals, can_fetch_across_multiple_contigs)
{
    // interval contig1:[3, 11) and contig2:[2, 7), test intervals on
    // different contigs, also test "tstart <=" and "tend >="
    DataSet ds{GenomicIntervalsTests::inputDir + "two_contigs.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 11}, {"contig2", 2, 7}};
    EXPECT_EQ(correct, result);
}

// Test various invalid AlignmentSets
TEST(BAM_DataSetGenomicIntervals, throws_on_missing_rname)
{
    // missing "rname"
    DataSet ds{GenomicIntervalsTests::inputDir + "invalid_missing_rname.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);
    EXPECT_THROW(ds.GenomicIntervals(), std::runtime_error);
}

TEST(BAM_DataSetGenomicIntervals, throws_on_invalid_rname_operator)
{
    // non-sensical "rname" operator ">"
    DataSet ds{GenomicIntervalsTests::inputDir + "invalid_rname_operator.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);
    EXPECT_THROW(ds.GenomicIntervals(), std::runtime_error);
}

TEST(BAM_DataSetGenomicIntervals, throws_on_invalid_tstart_operator)
{
    // non-sensical "tstart" operator "="
    DataSet ds{GenomicIntervalsTests::inputDir + "invalid_tstart_operator.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);
    EXPECT_THROW(ds.GenomicIntervals(), std::runtime_error);
}
