// Author: Derek Barnett

#include <cstddef>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/FastaWriter.h>
#include <pbbam/Unused.h>
#include <boost/algorithm/string.hpp>

#include "FastxTests.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace FastaReaderTests {

void CheckFastaSequence(const size_t index, const FastaSequence& seq)
{
    SCOPED_TRACE("checking FASTA seq:" + std::to_string(index));
    const auto& expected = FastxTests::ExpectedFasta.at(index);
    EXPECT_EQ(expected.Name(), seq.Name());
    EXPECT_EQ(expected.Bases(), seq.Bases());
}

}  // namespace FastaReaderTests

TEST(FastaReaderTest, throws_on_empty_filename)
{
    EXPECT_THROW(FastaReader reader{""}, std::runtime_error);
}

TEST(FastaReaderTest, throws_on_invalid_extension)
{
    EXPECT_THROW(FastaReader reader{"wrong.ext"}, std::runtime_error);
}

TEST(FastaReaderTest, can_open_text_fasta)
{
    const auto& fn = FastxTests::simpleFastaFn;
    EXPECT_NO_THROW(FastaReader reader{fn});
}

TEST(FastaReaderTest, can_open_gzip_fasta)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    EXPECT_NO_THROW(FastaReader reader{fn});
}

TEST(FastaReaderTest, can_open_bgzf_fasta)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    EXPECT_NO_THROW(FastaReader reader{fn});
}

TEST(FastaReaderTest, can_iterate_manually_on_text_fasta)
{
    const auto& fn = FastxTests::simpleFastaFn;
    FastaReader reader{fn};

    size_t count = 0;
    FastaSequence seq;
    while (reader.GetNext(seq)) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_iterate_manually_on_gzip_fasta)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    FastaReader reader{fn};

    size_t count = 0;
    FastaSequence seq;
    while (reader.GetNext(seq)) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_iterate_manually_on_bgzf_fasta)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    FastaReader reader{fn};

    size_t count = 0;
    FastaSequence seq;
    while (reader.GetNext(seq)) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_iterate_using_range_for_on_text_fasta)
{
    const auto& fn = FastxTests::simpleFastaFn;

    size_t count = 0;
    FastaReader reader{fn};
    for (const auto& seq : reader) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_iterate_using_range_for_on_gzip_fasta)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;

    size_t count = 0;
    FastaReader reader{fn};
    for (const auto& seq : reader) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_iterate_using_range_for_on_bgzf_fasta)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;

    size_t count = 0;
    FastaReader reader{fn};
    for (const auto& seq : reader) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_read_all_from_text_fasta)
{
    const auto& fn = FastxTests::simpleFastaFn;

    size_t count = 0;
    for (const auto& seq : FastaReader::ReadAll(fn)) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_read_all_from_gzip_fasta)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;

    size_t count = 0;
    for (const auto& seq : FastaReader::ReadAll(fn)) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_read_all_from_bgzf_fasta)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;

    size_t count = 0;
    for (const auto& seq : FastaReader::ReadAll(fn)) {
        FastaReaderTests::CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

TEST(FastaReaderTest, can_handle_windows_style_newlines)
{
    const std::string fn = FastxTests::fastxDataDir + "/windows_formatted.fasta";

    {
        size_t count = 0;
        FastaReader reader{fn};
        FastaSequence seq;
        while (reader.GetNext(seq)) {
            ++count;
            bool endOK = (boost::algorithm::ends_with(seq.Name(), "5p") ||
                          boost::algorithm::ends_with(seq.Name(), "3p"));
            EXPECT_TRUE(endOK);
        }
        EXPECT_EQ(7, count);  // 7 primers in total
    }
}
