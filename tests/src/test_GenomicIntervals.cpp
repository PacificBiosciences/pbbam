// Author: David Seifert

#include <pbbam/DataSet.h>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace GenomicIntervalsTests {

const std::string inputDir{PbbamTestsConfig::Data_Dir + "/test_GenomicIntervals/"};

}  // namespace GenomicIntervalsTests

TEST(DataSetGenomicIntervalsTest, NoFilter)
{
    // vanilla AlignmentSet, no filters
    DataSet ds{GenomicIntervalsTests::inputDir + "no_filter.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 0, 20}, {"contig2", 0, 10}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, Empty)
{
    // interval contig1:[5, 5), i.e., empty, yet both offsets are within range
    DataSet ds{GenomicIntervalsTests::inputDir + "empty.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct;
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, OutOfRange)
{
    // interval contig1:[1000, 10000), i.e., empty, as the selected range
    // lies above the contig1 size of 20
    DataSet ds{GenomicIntervalsTests::inputDir + "out_of_range.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct;
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, SingleInterval)
{
    // interval contig1:[3, 10)
    DataSet ds{GenomicIntervalsTests::inputDir + "single_interval.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 10}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, WholeContig)
{
    // interval contig1:[0, 20), i.e., select the whole contig
    DataSet ds{GenomicIntervalsTests::inputDir + "whole_contig.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 0, 20}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, ContigNameOnly)
{
    // interval contig1, i.e., select the whole contig, without a range filter
    DataSet ds{GenomicIntervalsTests::inputDir + "contig_name_only.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 0, 20}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, SingleIntervalLessOrEqual)
{
    // interval contig1:[3, 11), test "tstart <=" relation
    DataSet ds{GenomicIntervalsTests::inputDir + "single_interval_start_lte.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 11}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, SingleIntervalGreaterOrEqual)
{
    // interval contig1:[2, 10), test "tend >=" relation
    DataSet ds{GenomicIntervalsTests::inputDir + "single_interval_end_gte.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 2, 10}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, DisjointIntervals)
{
    // interval contig1:[3, 7),[13, 17), test that disjoint intervals remain disjoint
    DataSet ds{GenomicIntervalsTests::inputDir + "disjoint_intervals.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 7}, {"contig1", 13, 17}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, AdjacentIntervals)
{
    // interval contig1:[3, 17), test that intervals [3, 10) and [10, 17)
    // get merged into a single overall interval
    DataSet ds{GenomicIntervalsTests::inputDir + "adjacent_intervals.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);

    const std::vector<GenomicInterval> result = ds.GenomicIntervals();
    const std::vector<GenomicInterval> correct{{"contig1", 3, 17}};
    EXPECT_EQ(correct, result);
}

TEST(DataSetGenomicIntervalsTest, TwoContigs)
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
TEST(DataSetGenomicIntervalsTest, InvalidMissingRname)
{
    // missing "rname"
    DataSet ds{GenomicIntervalsTests::inputDir + "invalid_missing_rname.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);
    EXPECT_THROW(ds.GenomicIntervals(), std::runtime_error);
}

TEST(DataSetGenomicIntervalsTest, InvalidRnameOperator)
{
    // non-sensical "rname" operator ">"
    DataSet ds{GenomicIntervalsTests::inputDir + "invalid_rname_operator.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);
    EXPECT_THROW(ds.GenomicIntervals(), std::runtime_error);
}

TEST(DataSetGenomicIntervalsTest, InvalidTstartOperator)
{
    // non-sensical "tstart" operator "="
    DataSet ds{GenomicIntervalsTests::inputDir + "invalid_tstart_operator.alignmentset.xml"};
    ds.Type(DataSet::ALIGNMENT);
    EXPECT_THROW(ds.GenomicIntervals(), std::runtime_error);
}
