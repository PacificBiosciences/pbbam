// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#include <chrono>
#include <string>

#include <gtest/gtest.h>

#ifdef PBBAM_TESTING
#define private public
#endif

#include <pbbam/BamRecordBuilder.h>
#include <pbbam/BamTagCodec.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace tests {

static
void CheckRawData(const BamRecordImpl& bam)
{
    // ensure raw data (lengths at least) matches API-facing data

    const uint32_t expectedNameLength  = bam.Name().size() + 1;
    const uint32_t expectedNumCigarOps = bam.CigarData().size();
    const int32_t  expectedSeqLength   = bam.Sequence().length();
    const size_t   expectedTagsLength  = BamTagCodec::Encode(bam.Tags()).size();

    //  Name        CIGAR         Sequence       Quals      Tags
    // l_qname + (n_cigar * 4) + (l_qseq+1)/2 + l_qseq + << TAGS >>

    const int expectedTotalDataLength = expectedNameLength +
                                        (expectedNumCigarOps * 4) +
                                        (expectedSeqLength+1)/2 +
                                         expectedSeqLength +
                                         expectedTagsLength;

    EXPECT_TRUE((bool)bam.d_);
    EXPECT_EQ(expectedNameLength,      bam.d_->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps,     bam.d_->core.n_cigar);
    EXPECT_EQ(expectedSeqLength,       bam.d_->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, bam.d_->l_data);
}

static
void CheckRawData(const BamRecord& bam)
{ CheckRawData(bam.impl_); }

} // namespace tests

TEST(BamRecordBuilderTest, DefaultValues)
{
    BamRecordBuilder builder;
    BamRecord bam = builder.Build();

    const auto rawData = bam.impl_.d_;
    ASSERT_TRUE((bool)rawData);

    // fixed-length (core) data
    EXPECT_EQ(0, rawData->core.tid);
    EXPECT_EQ(0, rawData->core.pos);
    EXPECT_EQ(0, rawData->core.bin);
    EXPECT_EQ(0, rawData->core.qual);
    EXPECT_EQ(1, rawData->core.l_qname);    // initialized w/ NULL-term
    EXPECT_EQ(0, rawData->core.flag);
    EXPECT_EQ(0, rawData->core.n_cigar);
    EXPECT_EQ(0, rawData->core.l_qseq);
    EXPECT_EQ(0, rawData->core.mtid);
    EXPECT_EQ(0, rawData->core.mpos);
    EXPECT_EQ(0, rawData->core.isize);

    // variable length data
    EXPECT_TRUE(rawData->data != nullptr);
    EXPECT_EQ(1, rawData->l_data);
    EXPECT_EQ((int)0x800, rawData->m_data);  // check this if we change or tune later

    // -------------------------------
    // check data via API calls
    // -------------------------------

    EXPECT_EQ(0, bam.impl_.Bin());
    EXPECT_EQ(0, bam.impl_.Flag());
    EXPECT_EQ(0, bam.impl_.InsertSize());
    EXPECT_EQ(0, bam.impl_.MapQuality());
    EXPECT_EQ(0, bam.impl_.MateReferenceId());
    EXPECT_EQ(0, bam.impl_.MatePosition());
    EXPECT_EQ(0, bam.impl_.Position());
    EXPECT_EQ(0, bam.impl_.ReferenceId());
    EXPECT_EQ(0, bam.impl_.Tags().size());

    EXPECT_FALSE(bam.impl_.IsDuplicate());
    EXPECT_FALSE(bam.impl_.IsFailedQC());
    EXPECT_FALSE(bam.impl_.IsFirstMate());
    EXPECT_TRUE(bam.impl_.IsMapped());
    EXPECT_TRUE(bam.impl_.IsMateMapped());
    EXPECT_FALSE(bam.impl_.IsMateReverseStrand());
    EXPECT_FALSE(bam.impl_.IsPaired());
    EXPECT_TRUE(bam.impl_.IsPrimaryAlignment());
    EXPECT_FALSE(bam.impl_.IsProperPair());
    EXPECT_FALSE(bam.impl_.IsReverseStrand());
    EXPECT_FALSE(bam.impl_.IsSecondMate());
    EXPECT_FALSE(bam.impl_.IsSupplementaryAlignment());

    const std::string emptyString = "";
    EXPECT_EQ(emptyString, bam.impl_.Name());
    EXPECT_EQ(emptyString, bam.impl_.CigarData().ToStdString());
    EXPECT_EQ(emptyString, bam.impl_.Sequence());
    EXPECT_EQ(emptyString, bam.impl_.Qualities().Fastq());
    tests::CheckRawData(bam);
}

TEST(BamRecordBuilderTest, CheckSetters)
{
    // should be 28 bytes, encoded
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);

    BamRecordBuilder builder;
    builder.Bin(42)
           .Flag(42)
           .InsertSize(42)
           .MapQuality(42)
           .MatePosition(42)
           .MateReferenceId(42)
           .Position(42)
           .ReferenceId(42)
           .Tags(tags);

    BamRecord bam = builder.Build();

    // -------------------------------
    // check raw data
    // -------------------------------

    const auto rawData = bam.impl_.d_;
    ASSERT_TRUE((bool)rawData);

    // fixed-length (core) data
    EXPECT_EQ(42, rawData->core.tid);
    EXPECT_EQ(42, rawData->core.pos);
    EXPECT_EQ(42, rawData->core.bin);
    EXPECT_EQ(42, rawData->core.qual);
    EXPECT_EQ(1,  rawData->core.l_qname);    // initialized w/ NULL-term
    EXPECT_EQ(42, rawData->core.flag);
    EXPECT_EQ(0,  rawData->core.n_cigar);
    EXPECT_EQ(0,  rawData->core.l_qseq);
    EXPECT_EQ(42, rawData->core.mtid);
    EXPECT_EQ(42, rawData->core.mpos);
    EXPECT_EQ(42, rawData->core.isize);

    // variable length data
    EXPECT_TRUE(rawData->data != nullptr);
    EXPECT_EQ(29, rawData->l_data);         // NULL-term qname + tags
    EXPECT_EQ((int)0x800, rawData->m_data); // check this if we change or tune later

    // -------------------------------
    // check data via API calls
    // -------------------------------

    EXPECT_EQ(42, bam.impl_.Bin());
    EXPECT_EQ(42, bam.impl_.Flag());
    EXPECT_EQ(42, bam.impl_.InsertSize());
    EXPECT_EQ(42, bam.impl_.MapQuality());
    EXPECT_EQ(42, bam.impl_.MateReferenceId());
    EXPECT_EQ(42, bam.impl_.MatePosition());
    EXPECT_EQ(42, bam.impl_.Position());
    EXPECT_EQ(42, bam.impl_.ReferenceId());

    const TagCollection& fetchedTags = bam.impl_.Tags();

    EXPECT_TRUE(fetchedTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags.at("CA").ToUInt8Array());
}

//#define SEQ_LENGTH  7000
//#define NUM_RECORDS 1000

//const std::string& TEST_SEQUENCE  = std::string(SEQ_LENGTH, 'G');
//const std::string& TEST_QUALITIES = std::string(SEQ_LENGTH, '=');
//const std::string& TEST_NAME      = std::string(SEQ_LENGTH, '/');
//const std::string& TEST_TAGDATA   = std::string(SEQ_LENGTH, '2');

//TEST(BamRecordBuilderTest, JustDoingSomeTimings_BamRecordBuilder)
//{

//    BamRecordBuilder builder;

//    TagCollection tags;
//    tags["aa"] = TEST_TAGDATA;
//    tags["bb"] = TEST_TAGDATA;
//    tags["cc"] = TEST_TAGDATA;
//    tags["dd"] = TEST_TAGDATA;
//    tags["ee"] = TEST_TAGDATA;
//    tags["ff"] = TEST_TAGDATA;

//    auto start = std::chrono::steady_clock::now();

//    BamRecord record;
//    for (size_t i = 0; i < NUM_RECORDS; ++i) {
//        builder.Sequence(TEST_SEQUENCE)
//               .Qualities(TEST_QUALITIES)
//               .Name(TEST_NAME)
//               .Tags(tags)
//               .BuildInPlace(record);
//    }
//    auto end = std::chrono::steady_clock::now();
//    (void)record;
//    auto diff = end - start;
//    std::cout << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
//}


//TEST(BamRecordBuilderTest, JustDoingSomeTimings_BamRecordOnly)
//{
//    TagCollection tags;
//    tags["aa"] = TEST_TAGDATA;
//    tags["bb"] = TEST_TAGDATA;
//    tags["cc"] = TEST_TAGDATA;
//    tags["dd"] = TEST_TAGDATA;
//    tags["ee"] = TEST_TAGDATA;
//    tags["ff"] = TEST_TAGDATA;

//    auto start = std::chrono::steady_clock::now();

//    BamRecord record;
//    for (size_t i = 0; i < NUM_RECORDS; ++i) {
//        record.SetSequenceAndQualities(TEST_SEQUENCE, TEST_QUALITIES);
//        record.Name(TEST_NAME);
//        record.Tags(tags);
//    }
//    auto end = std::chrono::steady_clock::now();
//    (void)record;
//    auto diff = end - start;
//    std::cout << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
//}

