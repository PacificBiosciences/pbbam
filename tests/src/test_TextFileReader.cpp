// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include <pbbam/TextFileReader.h>

#include "FastxTests.h"

using TextFileReader = PacBio::BAM::TextFileReader;

// clang-format off

TEST(TextFileReaderTest, throws_on_empty_filename)
{
    EXPECT_THROW(TextFileReader reader{""}, std::runtime_error);
}

TEST(TextFileReaderTest, can_open_plain_text)
{
    const auto& fn = FastxTests::simpleFastaFn;
    EXPECT_NO_THROW(TextFileReader reader{fn});
}

TEST(TextFileReaderTest, can_open_gzip_text)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    EXPECT_NO_THROW(TextFileReader reader{fn});
}

TEST(TextFileReaderTest, can_open_bgzf_text)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    EXPECT_NO_THROW(TextFileReader reader{fn});
}

TEST(TextFileReaderTest, can_iterate_manually_on_plain_text)
{
    const auto& fn = FastxTests::simpleFastaFn;
    TextFileReader reader{fn};

    size_t count = 0;
    std::string line;
    while (reader.GetNext(line))
        ++count;

    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, count); // FASTA header + seq
}

TEST(TextFileReaderTest, can_iterate_manually_on_gzip_text)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    TextFileReader reader{fn};

    size_t count = 0;
    std::string line;
    while (reader.GetNext(line))
        ++count;

    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, count); // FASTA header + seq
}

TEST(TextFileReaderTest, can_iterate_manually_on_bgzf_text)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    TextFileReader reader{fn};

    size_t count = 0;
    std::string line;
    while (reader.GetNext(line))
        ++count;

    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, count); // FASTA header + seq
}

TEST(TextFileReaderTest, can_iterate_using_range_for_on_plain_text)
{
    const auto& fn = FastxTests::simpleFastaFn;
    TextFileReader reader{fn};

    size_t count = 0;
    for (const auto& line : reader) {
        EXPECT_FALSE(line.empty());
        ++count;
    }

    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, count); // FASTA header + seq
}

TEST(TextFileReaderTest, can_iterate_using_range_for_on_gzip_text)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    TextFileReader reader{fn};

    size_t count = 0;
    for (const auto& line : reader) {
        EXPECT_FALSE(line.empty());
        ++count;
    }

    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, count); // FASTA header + seq
}

TEST(TextFileReaderTest, can_iterate_using_range_for_on_bgzf_text)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    TextFileReader reader{fn};

    size_t count = 0;
    for (const auto& line : reader) {
        EXPECT_FALSE(line.empty());
        ++count;
    }

    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, count); // FASTA header + seq
}

TEST(TextFileReaderTest, can_read_all_from_plain_text)
{
    const auto& fn = FastxTests::simpleFastaFn;
    const auto lines = TextFileReader::ReadAll(fn);
    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, lines.size()); // FASTA header + seq)
}

TEST(TextFileReaderTest, can_read_all_from_gzip_text)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    const auto lines = TextFileReader::ReadAll(fn);
    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, lines.size()); // FASTA header + seq)
}

TEST(TextFileReaderTest, can_read_all_from_bgzf_text)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    const auto lines = TextFileReader::ReadAll(fn);
    EXPECT_EQ(FastxTests::ExpectedFasta.size() * 2, lines.size()); // FASTA header + seq)
}

// clang-foramt on
