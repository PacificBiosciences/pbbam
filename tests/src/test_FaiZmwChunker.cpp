// Author: Derek Barnett

#include "../../src/FaiZmwChunker.h"

#include <vector>

#include <gtest/gtest.h>

#include "FastxTests.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(FaiZmwChunkerTest, empty_input_zmws_yields_no_chunks)
{
    FaiIndex index;
    FaiZmwChunker chunker{index, 5};
    EXPECT_EQ(0, chunker.NumChunks());
}

TEST(FaiZmwChunkerTest, throws_if_requested_num_chunks_is_zero)
{
    FaiIndex index;
    EXPECT_THROW(FaiZmwChunker(index, 0), std::runtime_error);
}

TEST(FaiZmwChunkerTest, standard_chunking)
{
    {
        FaiZmwChunker chunker{FastxTests::chunkingFastaFaiFn, 5};
        ASSERT_EQ(5, chunker.NumChunks());

        // 7-7-6-6-6

        EXPECT_EQ(7, chunker.Chunk(0).NumZmws);
        EXPECT_EQ(7, chunker.Chunk(1).NumZmws);
        EXPECT_EQ(6, chunker.Chunk(2).NumZmws);
        EXPECT_EQ(6, chunker.Chunk(3).NumZmws);
        EXPECT_EQ(6, chunker.Chunk(4).NumZmws);

        EXPECT_EQ(7, chunker.Chunk(0).NumRecords);
        EXPECT_EQ(7, chunker.Chunk(1).NumRecords);
        EXPECT_EQ(6, chunker.Chunk(2).NumRecords);
        EXPECT_EQ(6, chunker.Chunk(3).NumRecords);
        EXPECT_EQ(9, chunker.Chunk(4).NumRecords);  // 3 records share zmw

        EXPECT_EQ("seq/0", chunker.Chunk(0).FirstSeqName);
        EXPECT_EQ("seq/7", chunker.Chunk(1).FirstSeqName);
        EXPECT_EQ("seq/14", chunker.Chunk(2).FirstSeqName);
        EXPECT_EQ("seq/20", chunker.Chunk(3).FirstSeqName);
        EXPECT_EQ("seq/50", chunker.Chunk(4).FirstSeqName);

        EXPECT_EQ(7, chunker.Chunk(0).FirstSeqOffset);
        EXPECT_EQ(91, chunker.Chunk(1).FirstSeqOffset);
        EXPECT_EQ(180, chunker.Chunk(2).FirstSeqOffset);
        EXPECT_EQ(258, chunker.Chunk(3).FirstSeqOffset);
        EXPECT_EQ(336, chunker.Chunk(4).FirstSeqOffset);
    }
    {
        FaiZmwChunker chunker{FastxTests::chunkingFastaFaiFn, 3};
        ASSERT_EQ(3, chunker.NumChunks());

        // 11-11-10

        EXPECT_EQ(11, chunker.Chunk(0).NumZmws);
        EXPECT_EQ(11, chunker.Chunk(1).NumZmws);
        EXPECT_EQ(10, chunker.Chunk(2).NumZmws);

        EXPECT_EQ(11, chunker.Chunk(0).NumRecords);
        EXPECT_EQ(11, chunker.Chunk(1).NumRecords);
        EXPECT_EQ(13, chunker.Chunk(2).NumRecords);  // 3 records share zmw

        EXPECT_EQ("seq/0", chunker.Chunk(0).FirstSeqName);
        EXPECT_EQ("seq/11", chunker.Chunk(1).FirstSeqName);
        EXPECT_EQ("seq/30", chunker.Chunk(2).FirstSeqName);

        EXPECT_EQ(7, chunker.Chunk(0).FirstSeqOffset);
        EXPECT_EQ(141, chunker.Chunk(1).FirstSeqOffset);
        EXPECT_EQ(284, chunker.Chunk(2).FirstSeqOffset);
    }
}

TEST(FaiZmwChunkerTest, one_chunk_contains_all_records)
{
    FaiZmwChunker chunker{FastxTests::chunkingFastaFaiFn, 1};
    ASSERT_EQ(1, chunker.NumChunks());

    // 32

    EXPECT_EQ(32, chunker.Chunk(0).NumZmws);
    EXPECT_EQ(35, chunker.Chunk(0).NumRecords);
    EXPECT_EQ("seq/0", chunker.Chunk(0).FirstSeqName);
    EXPECT_EQ(7, chunker.Chunk(0).FirstSeqOffset);
}

TEST(FaiZmwChunkerTest, one_zmw_per_chunk_if_requested_count_is_larger_than_input)
{
    FaiZmwChunker chunker{FastxTests::chunkingFastaFaiFn, 50};
    ASSERT_EQ(32, chunker.NumChunks());
    // 32 unique ZMWs
}
