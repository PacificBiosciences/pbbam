// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#ifdef PBBAM_TESTING
#define private public
#endif

#include <gtest/gtest.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamTagCodec.h>
#include <pbbam/Tag.h>
#include <pbbam/TagCollection.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;

namespace tests {

struct Bam1Deleter
{
    void operator()(bam1_t* b) {
        if (b)
            bam_destroy1(b);
        b = nullptr;
    }
};

static
BamRecord CreateRecord(void)
{
    BamRecord bam;
    bam.Bin(42);
    bam.Flag(42);
    bam.InsertSize(42);
    bam.MapQuality(42);
    bam.MatePosition(42);
    bam.MateReferenceId(42);
    bam.Position(42);
    bam.ReferenceId(42);

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam.Tags(tags);

    return bam;
}

static
void CheckRawData(const BamRecord& bam)
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

    EXPECT_TRUE((bool)bam.RawData());
    EXPECT_EQ(expectedNameLength,      bam.RawData()->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps,     bam.RawData()->core.n_cigar);
    EXPECT_EQ(expectedSeqLength,       bam.RawData()->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, bam.RawData()->l_data);
}

} // namespace tests

TEST(BamRecordCoreTest, RawDataDefaultValues)
{
    std::shared_ptr<bam1_t> rawData(bam_init1(), tests::Bam1Deleter());
    ASSERT_TRUE((bool)rawData);

    // fixed-length (core) data
    EXPECT_EQ(0, rawData->core.tid);
    EXPECT_EQ(0, rawData->core.pos);
    EXPECT_EQ(0, rawData->core.bin);
    EXPECT_EQ(0, rawData->core.qual);
    EXPECT_EQ(0, rawData->core.l_qname);
    EXPECT_EQ(0, rawData->core.flag);
    EXPECT_EQ(0, rawData->core.n_cigar);
    EXPECT_EQ(0, rawData->core.l_qseq);
    EXPECT_EQ(0, rawData->core.mtid);
    EXPECT_EQ(0, rawData->core.mpos);
    EXPECT_EQ(0, rawData->core.isize);

    // variable length data
    EXPECT_EQ(0, rawData->data);
    EXPECT_EQ(0, rawData->l_data);
    EXPECT_EQ(0, rawData->m_data);
}

TEST(BamRecordCoreTest, DefaultValues)
{
    BamRecord bam;

    // -------------------------------
    // check raw data
    // -------------------------------

    const std::shared_ptr<bam1_t> rawData = bam.RawData();
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

    EXPECT_EQ(0, bam.Bin());
    EXPECT_EQ(0, bam.Flag());
    EXPECT_EQ(0, bam.InsertSize());
    EXPECT_EQ(0, bam.MapQuality());
    EXPECT_EQ(0, bam.MateReferenceId());
    EXPECT_EQ(0, bam.MatePosition());
    EXPECT_EQ(0, bam.Position());
    EXPECT_EQ(0, bam.ReferenceId());
    EXPECT_EQ(0, bam.Tags().size());

    EXPECT_FALSE(bam.IsDuplicate());
    EXPECT_FALSE(bam.IsFailedQC());
    EXPECT_FALSE(bam.IsFirstMate());
    EXPECT_TRUE(bam.IsMapped());
    EXPECT_TRUE(bam.IsMateMapped());
    EXPECT_FALSE(bam.IsMateReverseStrand());
    EXPECT_FALSE(bam.IsPaired());
    EXPECT_TRUE(bam.IsPrimaryAlignment());
    EXPECT_FALSE(bam.IsProperPair());
    EXPECT_FALSE(bam.IsReverseStrand());
    EXPECT_FALSE(bam.IsSecondMate());
    EXPECT_FALSE(bam.IsSupplementaryAlignment());

    const std::string emptyString = "";
    EXPECT_EQ(emptyString, bam.Name());
    EXPECT_EQ(emptyString, bam.CigarData().ToStdString());
    EXPECT_EQ(emptyString, bam.Sequence());
    EXPECT_EQ(emptyString, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordCoreTest, CoreSetters)
{
    BamRecord bam;
    bam.Bin(42);
    bam.Flag(42);
    bam.InsertSize(42);
    bam.MapQuality(42);
    bam.MatePosition(42);
    bam.MateReferenceId(42);
    bam.Position(42);
    bam.ReferenceId(42);

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam.Tags(tags); // (28 bytes encoded)

    // -------------------------------
    // check raw data
    // -------------------------------

    const std::shared_ptr<bam1_t> rawData = bam.RawData();
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

    EXPECT_EQ(42, bam.Bin());
    EXPECT_EQ(42, bam.Flag());
    EXPECT_EQ(42, bam.InsertSize());
    EXPECT_EQ(42, bam.MapQuality());
    EXPECT_EQ(42, bam.MateReferenceId());
    EXPECT_EQ(42, bam.MatePosition());
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.ReferenceId());

    const TagCollection& fetchedTags = bam.Tags();

    EXPECT_TRUE(fetchedTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags.at("CA").ToUInt8Array());
}

TEST(BamRecordCoreTest, DeepCopyFromRawData)
{
    // init raw data
    std::shared_ptr<bam1_t> rawData(bam_init1(), tests::Bam1Deleter());
    ASSERT_TRUE((bool)rawData);

    rawData->core.tid = 42;
    rawData->core.pos = 42;
    rawData->core.bin = 42;
    rawData->core.qual = 42;
    rawData->core.flag = 42;
    rawData->core.mtid = 42;
    rawData->core.mpos = 42;
    rawData->core.isize = 42;

    const int32_t x = 42;
    char valueBytes[sizeof x];
    std::copy(static_cast<const char*>(static_cast<const void*>(&x)),
              static_cast<const char*>(static_cast<const void*>(&x)) + sizeof x,
              valueBytes);
    bam_aux_append(rawData.get(), "XY", 'i', sizeof(x), (uint8_t*)&valueBytes[0]);

    EXPECT_EQ(42, rawData->core.tid);
    EXPECT_EQ(42, rawData->core.pos);
    EXPECT_EQ(42, rawData->core.bin);
    EXPECT_EQ(42, rawData->core.qual);
    EXPECT_EQ(0,  rawData->core.l_qname);
    EXPECT_EQ(42, rawData->core.flag);
    EXPECT_EQ(0,  rawData->core.n_cigar);
    EXPECT_EQ(0,  rawData->core.l_qseq);
    EXPECT_EQ(42, rawData->core.mtid);
    EXPECT_EQ(42, rawData->core.mpos);
    EXPECT_EQ(42, rawData->core.isize);
    const int32_t fetchedX = bam_aux2i( bam_aux_get(rawData.get(), "XY") );
    EXPECT_EQ(42, fetchedX);

    // static "ctor"
    BamRecord bam = BamRecord::FromRawData(rawData);

    // make sure raw data is still valid
    EXPECT_EQ(42, rawData->core.tid);
    EXPECT_EQ(42, rawData->core.pos);
    EXPECT_EQ(42, rawData->core.bin);
    EXPECT_EQ(42, rawData->core.qual);
    EXPECT_EQ(0,  rawData->core.l_qname);
    EXPECT_EQ(42, rawData->core.flag);
    EXPECT_EQ(0,  rawData->core.n_cigar);
    EXPECT_EQ(0,  rawData->core.l_qseq);
    EXPECT_EQ(42, rawData->core.mtid);
    EXPECT_EQ(42, rawData->core.mpos);
    EXPECT_EQ(42, rawData->core.isize);
    EXPECT_TRUE(rawData->data != nullptr);
    EXPECT_TRUE(0 != rawData->l_data);
    EXPECT_TRUE(0 != rawData->m_data);

    // check new record
    EXPECT_EQ(42, bam.Bin());
    EXPECT_EQ(42, bam.Flag());
    EXPECT_EQ(42, bam.InsertSize());
    EXPECT_EQ(42, bam.MapQuality());
    EXPECT_EQ(42, bam.MateReferenceId());
    EXPECT_EQ(42, bam.MatePosition());
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.ReferenceId());
    EXPECT_EQ(x,  bam.Tags()["XY"].ToInt32());

    EXPECT_TRUE(bam.RawData()->data != nullptr);
    EXPECT_TRUE(bam.RawData()->m_data >= (int)0x800); // check this if we change or tune later

    // tweak raw data, make sure we've done a deep copy (so BamRecord isn't changed)
    rawData->core.pos = 37;
    EXPECT_EQ(37, rawData->core.pos);
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.RawData()->core.pos);
}

TEST(BamRecordCoreTest, CopyAssignment)
{
    BamRecord bam1;
    bam1.Bin(42);
    bam1.Flag(42);
    bam1.InsertSize(42);
    bam1.MapQuality(42);
    bam1.MatePosition(42);
    bam1.MateReferenceId(42);
    bam1.Position(42);
    bam1.ReferenceId(42);

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam1.Tags(tags);

    BamRecord bam2;
    bam2 = bam1;

    EXPECT_EQ(42, bam1.Bin());
    EXPECT_EQ(42, bam1.Flag());
    EXPECT_EQ(42, bam1.InsertSize());
    EXPECT_EQ(42, bam1.MapQuality());
    EXPECT_EQ(42, bam1.MateReferenceId());
    EXPECT_EQ(42, bam1.MatePosition());
    EXPECT_EQ(42, bam1.Position());
    EXPECT_EQ(42, bam1.ReferenceId());

    const TagCollection& fetchedTags1 = bam1.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    EXPECT_EQ(42, bam2.Bin());
    EXPECT_EQ(42, bam2.Flag());
    EXPECT_EQ(42, bam2.InsertSize());
    EXPECT_EQ(42, bam2.MapQuality());
    EXPECT_EQ(42, bam2.MateReferenceId());
    EXPECT_EQ(42, bam2.MatePosition());
    EXPECT_EQ(42, bam2.Position());
    EXPECT_EQ(42, bam2.ReferenceId());

    const TagCollection& fetchedTags2 = bam2.Tags();
    EXPECT_TRUE(fetchedTags2.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags2.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags2.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags2.at("CA").ToUInt8Array());

    tests::CheckRawData(bam1);
    tests::CheckRawData(bam2);
}

TEST(BamRecordCoreTest, SelfAssignmentTolerated)
{
    BamRecord bam1;
    bam1.Bin(42);
    bam1.Flag(42);
    bam1.InsertSize(42);
    bam1.MapQuality(42);
    bam1.MatePosition(42);
    bam1.MateReferenceId(42);
    bam1.Position(42);
    bam1.ReferenceId(42);

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam1.Tags(tags);

    bam1 = bam1;

    EXPECT_EQ(42, bam1.Bin());
    EXPECT_EQ(42, bam1.Flag());
    EXPECT_EQ(42, bam1.InsertSize());
    EXPECT_EQ(42, bam1.MapQuality());
    EXPECT_EQ(42, bam1.MateReferenceId());
    EXPECT_EQ(42, bam1.MatePosition());
    EXPECT_EQ(42, bam1.Position());
    EXPECT_EQ(42, bam1.ReferenceId());

    const TagCollection& fetchedTags1 = bam1.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    tests::CheckRawData(bam1);
}

TEST(BamRecordCoreTest, CopyConstructor)
{
    BamRecord bam1;
    bam1.Bin(42);
    bam1.Flag(42);
    bam1.InsertSize(42);
    bam1.MapQuality(42);
    bam1.MatePosition(42);
    bam1.MateReferenceId(42);
    bam1.Position(42);
    bam1.ReferenceId(42);

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam1.Tags(tags);

    BamRecord bam2(bam1);

    EXPECT_EQ(42, bam1.Bin());
    EXPECT_EQ(42, bam1.Flag());
    EXPECT_EQ(42, bam1.InsertSize());
    EXPECT_EQ(42, bam1.MapQuality());
    EXPECT_EQ(42, bam1.MateReferenceId());
    EXPECT_EQ(42, bam1.MatePosition());
    EXPECT_EQ(42, bam1.Position());
    EXPECT_EQ(42, bam1.ReferenceId());

    const TagCollection& fetchedTags1 = bam1.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    EXPECT_EQ(42, bam2.Bin());
    EXPECT_EQ(42, bam2.Flag());
    EXPECT_EQ(42, bam2.InsertSize());
    EXPECT_EQ(42, bam2.MapQuality());
    EXPECT_EQ(42, bam2.MateReferenceId());
    EXPECT_EQ(42, bam2.MatePosition());
    EXPECT_EQ(42, bam2.Position());
    EXPECT_EQ(42, bam2.ReferenceId());

    const TagCollection& fetchedTags2 = bam2.Tags();
    EXPECT_TRUE(fetchedTags2.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags2.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags2.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags2.at("CA").ToUInt8Array());

    tests::CheckRawData(bam1);
    tests::CheckRawData(bam2);
}

TEST(BamRecordCoreTest, CreateRecord_InternalTest)
{
    BamRecord bam = tests::CreateRecord();

    EXPECT_EQ(42, bam.Bin());
    EXPECT_EQ(42, bam.Flag());
    EXPECT_EQ(42, bam.InsertSize());
    EXPECT_EQ(42, bam.MapQuality());
    EXPECT_EQ(42, bam.MateReferenceId());
    EXPECT_EQ(42, bam.MatePosition());
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.ReferenceId());

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam.Tags(tags);

    tests::CheckRawData(bam);
}

TEST(BamRecordCoreTest, MoveAssignment)
{
    BamRecord bam;
    bam = std::move(tests::CreateRecord());

    EXPECT_EQ(42, bam.Bin());
    EXPECT_EQ(42, bam.Flag());
    EXPECT_EQ(42, bam.InsertSize());
    EXPECT_EQ(42, bam.MapQuality());
    EXPECT_EQ(42, bam.MateReferenceId());
    EXPECT_EQ(42, bam.MatePosition());
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.ReferenceId());

    const TagCollection& fetchedTags1 = bam.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    tests::CheckRawData(bam);
}

TEST(BamRecordCoreTest, MoveConstructor)
{
    BamRecord bam(std::move(tests::CreateRecord()));

    EXPECT_EQ(42, bam.Bin());
    EXPECT_EQ(42, bam.Flag());
    EXPECT_EQ(42, bam.InsertSize());
    EXPECT_EQ(42, bam.MapQuality());
    EXPECT_EQ(42, bam.MateReferenceId());
    EXPECT_EQ(42, bam.MatePosition());
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.ReferenceId());

    const TagCollection& fetchedTags1 = bam.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    tests::CheckRawData(bam);
}

TEST(BamRecordCoreTest, AlignmentFlags)
{
    // same set of flags, different ways of getting there

    // raw number
    BamRecord bam1;
    bam1.Flag(1107);

    // enum values
    BamRecord bam2;
    bam2.Flag(BamRecord::DUPLICATE |
              BamRecord::MATE_1 |
              BamRecord::REVERSE_STRAND |
              BamRecord::PROPER_PAIR |
              BamRecord::PAIRED
             );

    // convenience calls
    BamRecord bam3;
    bam3.SetDuplicate(true);
    bam3.SetFirstMate(true);
    bam3.SetReverseStrand(true);
    bam3.SetMapped(true);
    bam3.SetMateMapped(true);
    bam3.SetPaired(true);
    bam3.SetProperPair(true);
    bam3.SetPrimaryAlignment(true);

    // make sure all are same
    EXPECT_EQ(1107, bam1.Flag());
    EXPECT_EQ(1107, bam2.Flag());
    EXPECT_EQ(1107, bam3.Flag());

    // check API calls
    EXPECT_TRUE(bam1.IsPaired());
    EXPECT_TRUE(bam1.IsProperPair());
    EXPECT_TRUE(bam1.IsMapped());
    EXPECT_TRUE(bam1.IsMateMapped());
    EXPECT_TRUE(bam1.IsReverseStrand());
    EXPECT_FALSE(bam1.IsMateReverseStrand());
    EXPECT_TRUE(bam1.IsFirstMate());
    EXPECT_FALSE(bam1.IsSecondMate());
    EXPECT_TRUE(bam1.IsPrimaryAlignment());
    EXPECT_FALSE(bam1.IsFailedQC());
    EXPECT_TRUE(bam1.IsDuplicate());
    EXPECT_FALSE(bam1.IsSupplementaryAlignment());
}
