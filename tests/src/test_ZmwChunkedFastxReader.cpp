// Author: Derek Barnett

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "pbbam/ZmwChunkedFastaReader.h"
#include "pbbam/ZmwChunkedFastqReader.h"

#include "FastxTests.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(ZmwChunkedFastaReader, standard_fasta_from_chunk)
{
    ZmwChunkedFastaReader reader{FastxTests::chunkingFastaFn, 5};

    {
        const std::vector<std::string> expectedNames{"seq/0", "seq/1", "seq/2", "seq/3",
                                                     "seq/4", "seq/5", "seq/6"};
        reader.Chunk(0);

        std::vector<std::string> names;
        for (const auto& seq : reader)
            names.push_back(seq.Name());
        EXPECT_EQ(expectedNames, names);
    }
    {
        const std::vector<std::string> expectedNames{"seq/14", "seq/15", "seq/16",
                                                     "seq/17", "seq/18", "seq/19"};
        reader.Chunk(2);

        std::vector<std::string> names;
        for (const auto& seq : reader)
            names.push_back(seq.Name());
        EXPECT_EQ(expectedNames, names);
    }
    {
        const std::vector<std::string> expectedNames{
            "seq/50",          "seq/100/0_100",      "seq/100/100_200",
            "seq/100/200_300", "seq/100/300_400",    "seq/110/ccs",
            "seq/120/ccs",     "seq/130/transcript", "seq/140/transcript"};
        reader.Chunk(4);

        std::vector<std::string> names;
        for (const auto& seq : reader)
            names.push_back(seq.Name());
        EXPECT_EQ(expectedNames, names);
    }
}

TEST(ZmwChunkedFastqReader, standard_fastq_from_chunk)
{
    ZmwChunkedFastqReader reader{FastxTests::chunkingFastqFn, 5};

    {
        const std::vector<std::string> expectedNames{"seq/0", "seq/1", "seq/2", "seq/3",
                                                     "seq/4", "seq/5", "seq/6"};
        reader.Chunk(0);

        std::vector<std::string> names;
        for (const auto& seq : reader)
            names.push_back(seq.Name());
        EXPECT_EQ(expectedNames, names);
    }
    {
        const std::vector<std::string> expectedNames{"seq/14", "seq/15", "seq/16",
                                                     "seq/17", "seq/18", "seq/19"};
        reader.Chunk(2);

        std::vector<std::string> names;
        for (const auto& seq : reader)
            names.push_back(seq.Name());
        EXPECT_EQ(expectedNames, names);
    }
    {
        const std::vector<std::string> expectedNames{
            "seq/50",          "seq/100/0_100",      "seq/100/100_200",
            "seq/100/200_300", "seq/100/300_400",    "seq/110/ccs",
            "seq/120/ccs",     "seq/130/transcript", "seq/140/transcript"};
        reader.Chunk(4);

        std::vector<std::string> names;
        for (const auto& seq : reader)
            names.push_back(seq.Name());
        EXPECT_EQ(expectedNames, names);
    }
}