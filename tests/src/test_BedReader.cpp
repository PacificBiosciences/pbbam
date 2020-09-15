// Author: Derek Barnett

#include <pbbam/bed/BedReader.h>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/GenomicInterval.h>

#include "PbbamTestData.h"

using BedReader = PacBio::BAM::BedReader;
using GenomicInterval = PacBio::BAM::GenomicInterval;

// clang-format off

namespace BedReaderTests {

const std::string BedFn = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/bed/test.bed";
const std::string GzipBedFn = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/bed/test.bed.gz";

const std::vector<GenomicInterval> ExpectedIntervals {
    {"chr1", 213941196, 213942363},
    {"chr1", 213942363, 213943530},
    {"chr1", 213943530, 213944697},
    {"chr2", 158364697, 158365864},
    {"chr2", 158365864, 158367031},
    {"chr3", 127477031, 127478198},
    {"chr3", 127478198, 127479365},
    {"chr3", 127479365, 127480532},
    {"chr3", 127480532, 127481699}
};

void CheckManualIteration(const std::string& fn)
{
    size_t count = 0;
    BedReader reader{fn};
    GenomicInterval interval;
    while (reader.GetNext(interval)) {
        EXPECT_EQ(BedReaderTests::ExpectedIntervals.at(count), interval);
        ++count;
    }
    EXPECT_EQ(BedReaderTests::ExpectedIntervals.size(), count);
}

void CheckRangeFor(const std::string& fn)
{
    size_t count = 0;
    BedReader reader{fn};
    for (const auto& interval : reader) {
        EXPECT_EQ(BedReaderTests::ExpectedIntervals.at(count), interval);
        ++count;
    }
    EXPECT_EQ(BedReaderTests::ExpectedIntervals.size(), count);
}

void CheckReadAll(const std::string& fn)
{
    size_t count = 0;
    for (const auto& interval : BedReader::ReadAll(fn)) {
        EXPECT_EQ(BedReaderTests::ExpectedIntervals.at(count), interval);
        ++count;
    }
    EXPECT_EQ(BedReaderTests::ExpectedIntervals.size(), count);
}

}  // namespace BedReaderTests

TEST(BAM_BedReader, throws_on_empty_filename)
{
    EXPECT_THROW(BedReader reader{""}, std::runtime_error);
}

TEST(BAM_BedReader, throws_on_invalid_extension)
{
    EXPECT_THROW(BedReader reader{"wrong.ext"}, std::runtime_error);
}

TEST(BAM_BedReader, can_iterate_manually_on_text_bed)
{
    BedReaderTests::CheckManualIteration(BedReaderTests::BedFn);
}

TEST(BAM_BedReader, can_iterate_manually_on_gzip_bed)
{
    BedReaderTests::CheckManualIteration(BedReaderTests::GzipBedFn);
}

TEST(BAM_BedReader, can_iterate_using_range_for_on_text_bed)
{
    BedReaderTests::CheckRangeFor(BedReaderTests::BedFn);
}

TEST(BAM_BedReader, can_iterate_using_range_for_on_gzip_bed)
{
    BedReaderTests::CheckRangeFor(BedReaderTests::GzipBedFn);
}

TEST(BAM_BedReader, BedReaderTest_can_read_all_from_text_bed)
{
    BedReaderTests::CheckReadAll(BedReaderTests::BedFn);
}

TEST(BAM_BedReader, BedReaderTest_can_read_all_from_gzip_bed)
{
    BedReaderTests::CheckReadAll(BedReaderTests::GzipBedFn);
}

// clang-foramt on
