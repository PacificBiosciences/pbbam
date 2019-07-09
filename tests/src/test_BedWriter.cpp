// Author: Derek Barnett

#include <cstdio>

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/FormatUtils.h>
#include <pbbam/GenomicInterval.h>
#include <pbbam/bed/BedReader.h>
#include <pbbam/bed/BedWriter.h>

#include "PbbamTestData.h"

using BedReader = PacBio::BAM::BedReader;
using BedWriter = PacBio::BAM::BedWriter;
using GenomicInterval = PacBio::BAM::GenomicInterval;

namespace BedWriterTests {

const std::vector<GenomicInterval> Intervals{
    {"chr1", 213941196, 213942363}, {"chr1", 213942363, 213943530}, {"chr1", 213943530, 213944697},
    {"chr2", 158364697, 158365864}, {"chr2", 158365864, 158367031}, {"chr3", 127477031, 127478198},
    {"chr3", 127478198, 127479365}, {"chr3", 127479365, 127480532}, {"chr3", 127480532, 127481699}};

}  // namespace BedWriterTests

TEST(BedWriterTest, throws_on_empty_filename)
{
    EXPECT_THROW(BedWriter writer{""}, std::runtime_error);
}

TEST(BedWriterTest, can_write_plain_text)
{
    const std::string outFn = PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/out.bed";

    {
        BedWriter writer{outFn};
        for (const auto& interval : BedWriterTests::Intervals)
            writer.Write(interval);
    }
    EXPECT_EQ(PacBio::BAM::HtslibCompression::NONE,
              PacBio::BAM::FormatUtils::CompressionType(outFn));

    const auto contents = BedReader::ReadAll(outFn);
    EXPECT_TRUE(std::equal(BedWriterTests::Intervals.cbegin(), BedWriterTests::Intervals.cend(),
                           contents.cbegin()));

    remove(outFn.c_str());
}

TEST(BedWriterTest, can_write_gzipped_text)
{
    const std::string outFn = PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/out.bed.gz";

    {
        BedWriter writer{outFn};
        for (const auto& interval : BedWriterTests::Intervals)
            writer.Write(interval);
    }
    EXPECT_EQ(PacBio::BAM::HtslibCompression::GZIP,
              PacBio::BAM::FormatUtils::CompressionType(outFn));

    const auto contents = BedReader::ReadAll(outFn);
    EXPECT_TRUE(std::equal(BedWriterTests::Intervals.cbegin(), BedWriterTests::Intervals.cend(),
                           contents.cbegin()));

    remove(outFn.c_str());
}