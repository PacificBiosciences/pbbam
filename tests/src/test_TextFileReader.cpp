// Author: Derek Barnett

#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include <pbbam/TextFileReader.h>

#include "FastxTests.h"

using TextFileReader = PacBio::BAM::TextFileReader;

namespace TextFileReaderTests {

void CheckManualIteration(const std::string& fn)
{
    // FASTA header + seq
    const size_t expectedNumLines = FastxTests::ExpectedFasta.size() * 2;

    size_t count = 0;
    std::string line;
    TextFileReader reader{fn};
    while (reader.GetNext(line))
        ++count;
    EXPECT_EQ(expectedNumLines, count);
}

void CheckRangeFor(const std::string& fn)
{
    // FASTA header + seq
    const size_t expectedNumLines = FastxTests::ExpectedFasta.size() * 2;

    size_t count = 0;
    TextFileReader reader{fn};
    for (const auto& line : reader) {
        std::ignore = line;
        ++count;
    }
    EXPECT_EQ(expectedNumLines, count);
}

void CheckReadAll(const std::string& fn)
{
    // FASTA header + seq
    const size_t expectedNumLines = FastxTests::ExpectedFasta.size() * 2;

    size_t count = 0;
    for (const auto& line : TextFileReader::ReadAll(fn)) {
        std::ignore = line;
        ++count;
    }
    EXPECT_EQ(expectedNumLines, count);
}

}  // namespace TextFileReaderTests

// clang-format off

TEST(TextFileReaderTest, throws_on_empty_filename)
{
    EXPECT_THROW(TextFileReader reader{""}, std::runtime_error);
}

TEST(TextFileReaderTest, can_open_plain_text)
{
    EXPECT_NO_THROW(TextFileReader reader{FastxTests::simpleFastaFn});
}

TEST(TextFileReaderTest, can_open_gzip_text)
{
    EXPECT_NO_THROW(TextFileReader reader{FastxTests::simpleFastaGzipFn});
}

TEST(TextFileReaderTest, can_open_bgzf_text)
{
    EXPECT_NO_THROW(TextFileReader reader{FastxTests::simpleFastaBgzfFn});
}

TEST(TextFileReaderTest, can_iterate_manually_on_plain_text)
{
    TextFileReaderTests::CheckManualIteration(FastxTests::simpleFastaFn);
}

TEST(TextFileReaderTest, can_iterate_manually_on_gzip_text)
{
    TextFileReaderTests::CheckManualIteration(FastxTests::simpleFastaGzipFn);
}

TEST(TextFileReaderTest, can_iterate_manually_on_bgzf_text)
{
    TextFileReaderTests::CheckManualIteration(FastxTests::simpleFastaBgzfFn);
}

TEST(TextFileReaderTest, can_iterate_using_range_for_on_plain_text)
{
    TextFileReaderTests::CheckRangeFor(FastxTests::simpleFastaFn);
}

TEST(TextFileReaderTest, can_iterate_using_range_for_on_gzip_text)
{
    TextFileReaderTests::CheckRangeFor(FastxTests::simpleFastaGzipFn);
}

TEST(TextFileReaderTest, can_iterate_using_range_for_on_bgzf_text)
{
    TextFileReaderTests::CheckRangeFor(FastxTests::simpleFastaBgzfFn);
}

TEST(TextFileReaderTest, can_read_all_from_plain_text)
{
    TextFileReaderTests::CheckReadAll(FastxTests::simpleFastaFn);
}

TEST(TextFileReaderTest, can_read_all_from_gzip_text)
{
    TextFileReaderTests::CheckReadAll(FastxTests::simpleFastaGzipFn);
}

TEST(TextFileReaderTest, can_read_all_from_bgzf_text)
{
    TextFileReaderTests::CheckReadAll(FastxTests::simpleFastaBgzfFn);
}

// clang-foramt on
