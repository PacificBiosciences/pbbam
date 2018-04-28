// Author: Derek Barnett

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/BamRecordBuilder.h>
#include <pbbam/BamTagCodec.h>
#include "../src/MemoryUtils.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace BamRecordBuilderTests {

static void CheckRawData(const BamRecordImpl& bam)
{
    // ensure raw data (lengths at least) matches API-facing data

    const uint32_t expectedNameLength = bam.Name().size() + 1;
    const uint32_t expectedNumCigarOps = bam.CigarData().size();
    const int32_t expectedSeqLength = bam.Sequence().length();
    const size_t expectedTagsLength = BamTagCodec::Encode(bam.Tags()).size();

    //  Name        CIGAR         Sequence       Quals      Tags
    // l_qname + (n_cigar * 4) + (l_qseq+1)/2 + l_qseq + << TAGS >>

    const int expectedTotalDataLength = expectedNameLength + (expectedNumCigarOps * 4) +
                                        (expectedSeqLength + 1) / 2 + expectedSeqLength +
                                        expectedTagsLength;

    const auto rawData = PacBio::BAM::internal::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));
    EXPECT_EQ(expectedNameLength, rawData->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps, rawData->core.n_cigar);
    EXPECT_EQ(expectedSeqLength, rawData->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, rawData->l_data);
}

static void CheckRawData(const BamRecord& bam) { CheckRawData(bam.Impl()); }

}  // namespace BamRecordBuilderTests

TEST(BamRecordBuilderTest, DefaultValues)
{
    BamRecordBuilder builder;
    BamRecord bam = builder.Build();

    const auto rawData = PacBio::BAM::internal::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));

    // fixed-length (core) data
    EXPECT_EQ(0, rawData->core.tid);
    EXPECT_EQ(0, rawData->core.pos);
    EXPECT_EQ(0, rawData->core.bin);
    EXPECT_EQ(0, rawData->core.qual);
    EXPECT_EQ(1, rawData->core.l_qname);  // initialized w/ NULL-term
    EXPECT_EQ(0, rawData->core.flag);
    EXPECT_EQ(0, rawData->core.n_cigar);
    EXPECT_EQ(0, rawData->core.l_qseq);
    EXPECT_EQ(0, rawData->core.mtid);
    EXPECT_EQ(0, rawData->core.mpos);
    EXPECT_EQ(0, rawData->core.isize);

    // variable length data
    EXPECT_TRUE(rawData->data != nullptr);
    EXPECT_EQ(1, rawData->l_data);
    EXPECT_EQ(int{0x800}, rawData->m_data);  // check this if we change or tune later

    // -------------------------------
    // check data via API calls
    // -------------------------------

    EXPECT_EQ(0, bam.Impl().Bin());
    EXPECT_EQ(0, bam.Impl().Flag());
    EXPECT_EQ(0, bam.Impl().InsertSize());
    EXPECT_EQ(0, bam.Impl().MapQuality());
    EXPECT_EQ(0, bam.Impl().MateReferenceId());
    EXPECT_EQ(0, bam.Impl().MatePosition());
    EXPECT_EQ(0, bam.Impl().Position());
    EXPECT_EQ(0, bam.Impl().ReferenceId());
    EXPECT_EQ(0, bam.Impl().Tags().size());

    EXPECT_FALSE(bam.Impl().IsDuplicate());
    EXPECT_FALSE(bam.Impl().IsFailedQC());
    EXPECT_FALSE(bam.Impl().IsFirstMate());
    EXPECT_TRUE(bam.Impl().IsMapped());
    EXPECT_TRUE(bam.Impl().IsMateMapped());
    EXPECT_FALSE(bam.Impl().IsMateReverseStrand());
    EXPECT_FALSE(bam.Impl().IsPaired());
    EXPECT_TRUE(bam.Impl().IsPrimaryAlignment());
    EXPECT_FALSE(bam.Impl().IsProperPair());
    EXPECT_FALSE(bam.Impl().IsReverseStrand());
    EXPECT_FALSE(bam.Impl().IsSecondMate());
    EXPECT_FALSE(bam.Impl().IsSupplementaryAlignment());

    const std::string emptyString = "";
    EXPECT_EQ(emptyString, bam.Impl().Name());
    EXPECT_EQ(emptyString, bam.Impl().CigarData().ToStdString());
    EXPECT_EQ(emptyString, bam.Impl().Sequence());
    EXPECT_EQ(emptyString, bam.Impl().Qualities().Fastq());
    BamRecordBuilderTests::CheckRawData(bam);
}

TEST(BamRecordBuilderTest, CheckSetters)
{
    // should be 28 bytes, encoded
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

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

    const auto rawData = PacBio::BAM::internal::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));

    // fixed-length (core) data
    EXPECT_EQ(42, rawData->core.tid);
    EXPECT_EQ(42, rawData->core.pos);
    EXPECT_EQ(42, rawData->core.bin);
    EXPECT_EQ(42, rawData->core.qual);
    EXPECT_EQ(1, rawData->core.l_qname);  // initialized w/ NULL-term
    EXPECT_EQ(42, rawData->core.flag);
    EXPECT_EQ(0, rawData->core.n_cigar);
    EXPECT_EQ(0, rawData->core.l_qseq);
    EXPECT_EQ(42, rawData->core.mtid);
    EXPECT_EQ(42, rawData->core.mpos);
    EXPECT_EQ(42, rawData->core.isize);

    // variable length data
    EXPECT_TRUE(rawData->data != nullptr);
    EXPECT_EQ(29, rawData->l_data);          // NULL-term qname + tags
    EXPECT_EQ(int{0x800}, rawData->m_data);  // check this if we change or tune later

    // -------------------------------
    // check data via API calls
    // -------------------------------

    EXPECT_EQ(42, bam.Impl().Bin());
    EXPECT_EQ(42, bam.Impl().Flag());
    EXPECT_EQ(42, bam.Impl().InsertSize());
    EXPECT_EQ(42, bam.Impl().MapQuality());
    EXPECT_EQ(42, bam.Impl().MateReferenceId());
    EXPECT_EQ(42, bam.Impl().MatePosition());
    EXPECT_EQ(42, bam.Impl().Position());
    EXPECT_EQ(42, bam.Impl().ReferenceId());

    const TagCollection fetchedTags = bam.Impl().Tags();

    EXPECT_TRUE(fetchedTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags.at("CA").ToUInt8Array());
}
