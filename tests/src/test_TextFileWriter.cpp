#include <pbbam/TextFileWriter.h>

#include <cstdio>

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/FormatUtils.h>
#include <pbbam/TextFileReader.h>

#include "PbbamTestData.h"

using TextFileReader = PacBio::BAM::TextFileReader;
using TextFileWriter = PacBio::BAM::TextFileWriter;

namespace TextFileWriterTests {

const std::vector<std::string> Lines{"foo", "bar", "baz"};

void CheckRoundTrip(const std::string& outFn, const PacBio::BAM::HtslibCompression compressionType)
{
    {
        TextFileWriter writer{outFn};
        for (const auto& interval : TextFileWriterTests::Lines)
            writer.Write(interval);
    }
    EXPECT_EQ(compressionType, PacBio::BAM::FormatUtils::CompressionType(outFn));

    const auto contents = TextFileReader::ReadAll(outFn);
    EXPECT_TRUE(std::equal(TextFileWriterTests::Lines.cbegin(), TextFileWriterTests::Lines.cend(),
                           contents.cbegin()));

    remove(outFn.c_str());
}

}  // namespace TextFileWriterTests

TEST(BAM_TextFileWriter, throws_on_empty_filename)
{
    EXPECT_THROW(TextFileWriter writer{""}, std::runtime_error);
}

TEST(BAM_TextFileWriter, can_write_plain_text)
{
    TextFileWriterTests::CheckRoundTrip(
        PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/out.txt",
        PacBio::BAM::HtslibCompression::NONE);
}

TEST(BAM_TextFileWriter, can_write_gzipped_text)
{
    TextFileWriterTests::CheckRoundTrip(
        PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/out.txt.gz",
        PacBio::BAM::HtslibCompression::GZIP);
}
