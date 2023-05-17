#include <pbbam/FastaReader.h>

#include <cstddef>

#include <gtest/gtest.h>

#include <boost/algorithm/string.hpp>

#include "PbbamTestData.h"

#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/FastaWriter.h>

#include "FastxTests.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace FastaReaderTests {

void CheckFastaSequence(const std::size_t index, const FastaSequence& seq)
{
    SCOPED_TRACE("checking FASTA seq:" + std::to_string(index));
    const auto& expected = FastxTests::ExpectedFasta.at(index);
    EXPECT_EQ(expected.Name(), seq.Name());
    EXPECT_EQ(expected.Bases(), seq.Bases());
}

void CheckManualIteration(const std::string& fn)
{
    FastaReader reader{fn};

    std::size_t count = 0;
    FastaSequence seq;
    while (reader.GetNext(seq)) {
        CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

void CheckRangeFor(const std::string& fn)
{
    std::size_t count = 0;
    FastaReader reader{fn};
    for (const auto& seq : reader) {
        CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

void CheckReadAll(const std::string& fn)
{
    std::size_t count = 0;
    for (const auto& seq : FastaReader::ReadAll(fn)) {
        CheckFastaSequence(count, seq);
        ++count;
    }
    EXPECT_EQ(FastxTests::ExpectedFasta.size(), count);
}

}  // namespace FastaReaderTests

TEST(BAM_FastaReader, throws_on_empty_filename)
{
    EXPECT_THROW(FastaReader reader{""}, std::runtime_error);
}

TEST(BAM_FastaReader, throws_on_invalid_extension)
{
    EXPECT_THROW(FastaReader reader{"wrong.ext"}, std::runtime_error);
}

TEST(BAM_FastaReader, can_open_text_fasta)
{
    const auto& fn = FastxTests::simpleFastaFn;
    EXPECT_NO_THROW(FastaReader reader{fn});
}

TEST(BAM_FastaReader, can_open_gzip_fasta)
{
    const auto& fn = FastxTests::simpleFastaGzipFn;
    EXPECT_NO_THROW(FastaReader reader{fn});
}

TEST(BAM_FastaReader, can_open_bgzf_fasta)
{
    const auto& fn = FastxTests::simpleFastaBgzfFn;
    EXPECT_NO_THROW(FastaReader reader{fn});
}

TEST(BAM_FastaReader, can_iterate_manually_on_text_fasta)
{
    FastaReaderTests::CheckManualIteration(FastxTests::simpleFastaFn);
}

TEST(BAM_FastaReader, can_iterate_manually_on_text_fsa)
{
    FastaReaderTests::CheckManualIteration(FastxTests::simpleFsaFn);
}

TEST(BAM_FastaReader, can_iterate_manually_on_gzip_fasta)
{
    FastaReaderTests::CheckManualIteration(FastxTests::simpleFastaGzipFn);
}

TEST(BAM_FastaReader, can_iterate_manually_on_bgzf_fasta)
{
    FastaReaderTests::CheckManualIteration(FastxTests::simpleFastaBgzfFn);
}

TEST(BAM_FastaReader, can_iterate_using_range_for_on_text_fasta)
{
    FastaReaderTests::CheckRangeFor(FastxTests::simpleFastaFn);
}

TEST(BAM_FastaReader, can_iterate_using_range_for_on_gzip_fasta)
{
    FastaReaderTests::CheckRangeFor(FastxTests::simpleFastaGzipFn);
}

TEST(BAM_FastaReader, can_iterate_using_range_for_on_bgzf_fasta)
{
    FastaReaderTests::CheckRangeFor(FastxTests::simpleFastaBgzfFn);
}

TEST(BAM_FastaReader, can_read_all_from_text_fasta)
{
    FastaReaderTests::CheckReadAll(FastxTests::simpleFastaFn);
}

TEST(BAM_FastaReader, can_read_all_from_gzip_fasta)
{
    FastaReaderTests::CheckReadAll(FastxTests::simpleFastaGzipFn);
}

TEST(BAM_FastaReader, can_read_all_from_bgzf_fasta)
{
    FastaReaderTests::CheckReadAll(FastxTests::simpleFastaBgzfFn);
}

TEST(BAM_FastaReader, can_handle_windows_style_newlines)
{
    std::size_t count = 0;
    FastaReader reader{FastxTests::fastxDataDir + "/windows_formatted.fasta"};
    FastaSequence seq;
    while (reader.GetNext(seq)) {
        ++count;
        const bool endOK = (boost::algorithm::ends_with(seq.Name(), "5p") ||
                            boost::algorithm::ends_with(seq.Name(), "3p"));
        EXPECT_TRUE(endOK);
    }
    EXPECT_EQ(7, count);  // 7 primers in total
}
