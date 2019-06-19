// Author: Derek Barnett

#include <cstdio>

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/FormatUtils.h>
#include <pbbam/TextFileReader.h>
#include <pbbam/TextFileWriter.h>

#include "PbbamTestData.h"

using TextFileReader = PacBio::BAM::TextFileReader;
using TextFileWriter = PacBio::BAM::TextFileWriter;

TEST(TextFileWriterTest, throws_on_empty_filename)
{
    EXPECT_THROW(TextFileWriter writer{""}, std::runtime_error);
}

TEST(TextFileWriterTest, can_write_plain_text)
{
    const std::string outFn = PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/out.txt";
    const std::vector<std::string> lines{"foo", "bar", "baz"};

    {
        TextFileWriter writer{outFn};
        for (const auto& line : lines)
            writer.Write(line);
    }
    EXPECT_EQ(PacBio::BAM::HtslibCompression::NONE,
              PacBio::BAM::FormatUtils::CompressionType(outFn));

    const auto contents = TextFileReader::ReadAll(outFn);
    EXPECT_TRUE(std::equal(lines.cbegin(), lines.cend(), contents.cbegin()));

    remove(outFn.c_str());
}

TEST(TextFileWriterTest, can_write_gzipped_text)
{
    const std::string outFn = PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/out.txt.gz";
    const std::vector<std::string> lines{"foo", "bar", "baz"};

    {
        TextFileWriter writer{outFn};
        for (const auto& line : lines)
            writer.Write(line);
    }
    EXPECT_EQ(PacBio::BAM::HtslibCompression::GZIP,
              PacBio::BAM::FormatUtils::CompressionType(outFn));

    const auto contents = TextFileReader::ReadAll(outFn);
    EXPECT_TRUE(std::equal(lines.cbegin(), lines.cend(), contents.cbegin()));

    remove(outFn.c_str());
}
