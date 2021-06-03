#include <pbbam/BamRecordImpl.h>

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamTagCodec.h>
#include <pbbam/Tag.h>
#include <pbbam/TagCollection.h>

#include "../src/MemoryUtils.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace BamRecordImplCoreTests {

struct Bam1Deleter
{
    void operator()(bam1_t* b) const noexcept
    {
        if (b != nullptr) {
            bam_destroy1(b);
        }
        b = nullptr;
    }
};

static BamRecordImpl CreateBamImpl()
{
    BamRecordImpl bam;
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
    tags["XY"] = int32_t{-42};
    bam.Tags(tags);

    return bam;
}

static void CheckRawData(const BamRecordImpl& bam)
{
    // ensure raw data (lengths at least) matches API-facing data
    const uint32_t expectedNameBytes = bam.Name().size() + 1;  // include NULL term
    const uint32_t expectedNameNulls = 4 - (expectedNameBytes % 4);
    const uint32_t expectedNameLength = expectedNameBytes + expectedNameNulls;
    const uint32_t expectedNumCigarOps = bam.CigarData().size();
    const int32_t expectedSeqLength = bam.Sequence().length();
    const size_t expectedTagsLength = BamTagCodec::Encode(bam.Tags()).size();

    //  Name        CIGAR         Sequence       Quals      Tags
    // l_qname + (n_cigar * 4) + (l_qseq+1)/2 + l_qseq + << TAGS >>
    const int expectedTotalDataLength = expectedNameLength + (expectedNumCigarOps * 4) +
                                        (expectedSeqLength + 1) / 2 + expectedSeqLength +
                                        expectedTagsLength;

    const auto& rawData = PacBio::BAM::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));

    EXPECT_EQ(expectedNameNulls, rawData->core.l_extranul);
    EXPECT_EQ(expectedNameLength, rawData->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps, rawData->core.n_cigar);
    EXPECT_EQ(expectedSeqLength, rawData->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, rawData->l_data);
}

}  // namespace BamRecordImplCoreTests

TEST(BAM_BamRecordImplCore, initialized_with_correct_raw_htslib_values)
{
    std::shared_ptr<bam1_t> rawData(bam_init1(), BamRecordImplCoreTests::Bam1Deleter());
    ASSERT_TRUE(static_cast<bool>(rawData));

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
    EXPECT_EQ(0, rawData->l_data);  // initial aligned QNAME
    EXPECT_EQ(0, rawData->m_data);  // check this if we change or tune later
}

TEST(BAM_BamRecordImplCore, initialized_with_correct_pbbam_values)
{
    BamRecordImpl bam;

    // -------------------------------
    // check raw data
    // -------------------------------

    const auto& rawData = PacBio::BAM::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));

    // fixed-length (core) data
    // (forced init unmapped, with NULL-term as QNAME)
    EXPECT_EQ(-1, rawData->core.tid);
    EXPECT_EQ(-1, rawData->core.pos);
    EXPECT_EQ(0, rawData->core.bin);
    EXPECT_EQ(255, rawData->core.qual);
    EXPECT_EQ(3, rawData->core.l_extranul);  // alignment nulls
    EXPECT_EQ(4, rawData->core.l_qname);     // normal null term + alignment nulls
    EXPECT_EQ(BamRecordImpl::UNMAPPED, rawData->core.flag);
    EXPECT_EQ(0, rawData->core.n_cigar);
    EXPECT_EQ(0, rawData->core.l_qseq);
    EXPECT_EQ(-1, rawData->core.mtid);
    EXPECT_EQ(-1, rawData->core.mpos);
    EXPECT_EQ(0, rawData->core.isize);

    // variable length data
    EXPECT_TRUE(rawData->data != nullptr);
    EXPECT_EQ(4, rawData->l_data);  // initial aligned QNAME

    // -------------------------------
    // check data via API calls
    // -------------------------------

    EXPECT_EQ(0, bam.Bin());
    EXPECT_EQ(BamRecordImpl::UNMAPPED, bam.Flag());
    EXPECT_EQ(0, bam.InsertSize());
    EXPECT_EQ(255, bam.MapQuality());
    EXPECT_EQ(-1, bam.MateReferenceId());
    EXPECT_EQ(-1, bam.MatePosition());
    EXPECT_EQ(-1, bam.Position());
    EXPECT_EQ(-1, bam.ReferenceId());
    EXPECT_EQ(0, bam.Tags().size());

    EXPECT_FALSE(bam.IsDuplicate());
    EXPECT_FALSE(bam.IsFailedQC());
    EXPECT_FALSE(bam.IsFirstMate());
    EXPECT_FALSE(bam.IsMapped());
    EXPECT_TRUE(bam.IsMateMapped());
    EXPECT_FALSE(bam.IsMateReverseStrand());
    EXPECT_FALSE(bam.IsPaired());
    EXPECT_TRUE(bam.IsPrimaryAlignment());
    EXPECT_FALSE(bam.IsProperPair());
    EXPECT_FALSE(bam.IsReverseStrand());
    EXPECT_FALSE(bam.IsSecondMate());
    EXPECT_FALSE(bam.IsSupplementaryAlignment());

    const std::string emptyString;
    EXPECT_EQ(emptyString, bam.Name());
    EXPECT_EQ(emptyString, bam.CigarData().ToStdString());
    EXPECT_EQ(emptyString, bam.Sequence());
    EXPECT_EQ(emptyString, bam.Qualities().Fastq());
    BamRecordImplCoreTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplCore, can_be_modified_with_general_setters)
{
    BamRecordImpl bam;
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
    tags["XY"] = int32_t{-42};
    bam.Tags(tags);  // (28 bytes encoded)

    // -------------------------------
    // check raw data
    // -------------------------------

    const auto& rawData = PacBio::BAM::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));

    // fixed-length (core) data
    EXPECT_EQ(42, rawData->core.tid);
    EXPECT_EQ(42, rawData->core.pos);
    EXPECT_EQ(42, rawData->core.bin);
    EXPECT_EQ(42, rawData->core.qual);
    EXPECT_EQ(3, rawData->core.l_extranul);  // alignment nulls
    EXPECT_EQ(4, rawData->core.l_qname);     // normal null term + alignment nulls
    EXPECT_EQ(42, rawData->core.flag);
    EXPECT_EQ(0, rawData->core.n_cigar);
    EXPECT_EQ(0, rawData->core.l_qseq);
    EXPECT_EQ(42, rawData->core.mtid);
    EXPECT_EQ(42, rawData->core.mpos);
    EXPECT_EQ(42, rawData->core.isize);

    // variable length data
    EXPECT_TRUE(rawData->data != nullptr);
    EXPECT_EQ(32, rawData->l_data);  // aligned qname + tags

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

    const TagCollection fetchedTags = bam.Tags();

    EXPECT_TRUE(fetchedTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags.at("CA").ToUInt8Array());
}

TEST(BAM_BamRecordImplCore, can_be_copy_assigned)
{
    BamRecordImpl bam1;
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
    tags["XY"] = int32_t{-42};
    bam1.Tags(tags);

    BamRecordImpl bam2;
    bam2 = bam1;

    EXPECT_EQ(42, bam1.Bin());
    EXPECT_EQ(42, bam1.Flag());
    EXPECT_EQ(42, bam1.InsertSize());
    EXPECT_EQ(42, bam1.MapQuality());
    EXPECT_EQ(42, bam1.MateReferenceId());
    EXPECT_EQ(42, bam1.MatePosition());
    EXPECT_EQ(42, bam1.Position());
    EXPECT_EQ(42, bam1.ReferenceId());

    const TagCollection fetchedTags1 = bam1.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    EXPECT_EQ(42, bam2.Bin());
    EXPECT_EQ(42, bam2.Flag());
    EXPECT_EQ(42, bam2.InsertSize());
    EXPECT_EQ(42, bam2.MapQuality());
    EXPECT_EQ(42, bam2.MateReferenceId());
    EXPECT_EQ(42, bam2.MatePosition());
    EXPECT_EQ(42, bam2.Position());
    EXPECT_EQ(42, bam2.ReferenceId());

    const TagCollection fetchedTags2 = bam2.Tags();
    EXPECT_TRUE(fetchedTags2.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags2.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags2.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags2.at("CA").ToUInt8Array());

    BamRecordImplCoreTests::CheckRawData(bam1);
    BamRecordImplCoreTests::CheckRawData(bam2);
}

TEST(BAM_BamRecordImplCore, self_assignment_is_tolerated)
{
    BamRecordImpl bam1;
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
    tags["XY"] = int32_t{-42};
    bam1.Tags(tags);

    EXPECT_EQ(42, bam1.Bin());
    EXPECT_EQ(42, bam1.Flag());
    EXPECT_EQ(42, bam1.InsertSize());
    EXPECT_EQ(42, bam1.MapQuality());
    EXPECT_EQ(42, bam1.MateReferenceId());
    EXPECT_EQ(42, bam1.MatePosition());
    EXPECT_EQ(42, bam1.Position());
    EXPECT_EQ(42, bam1.ReferenceId());

    const TagCollection fetchedTags1 = bam1.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    BamRecordImplCoreTests::CheckRawData(bam1);
}

TEST(BAM_BamRecordImplCore, can_be_copy_constructed)
{
    BamRecordImpl bam1;
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
    tags["XY"] = int32_t{-42};
    bam1.Tags(tags);

    BamRecordImpl bam2(bam1);

    EXPECT_EQ(42, bam1.Bin());
    EXPECT_EQ(42, bam1.Flag());
    EXPECT_EQ(42, bam1.InsertSize());
    EXPECT_EQ(42, bam1.MapQuality());
    EXPECT_EQ(42, bam1.MateReferenceId());
    EXPECT_EQ(42, bam1.MatePosition());
    EXPECT_EQ(42, bam1.Position());
    EXPECT_EQ(42, bam1.ReferenceId());

    const TagCollection fetchedTags1 = bam1.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    EXPECT_EQ(42, bam2.Bin());
    EXPECT_EQ(42, bam2.Flag());
    EXPECT_EQ(42, bam2.InsertSize());
    EXPECT_EQ(42, bam2.MapQuality());
    EXPECT_EQ(42, bam2.MateReferenceId());
    EXPECT_EQ(42, bam2.MatePosition());
    EXPECT_EQ(42, bam2.Position());
    EXPECT_EQ(42, bam2.ReferenceId());

    const TagCollection fetchedTags2 = bam2.Tags();
    EXPECT_TRUE(fetchedTags2.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags2.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags2.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags2.at("CA").ToUInt8Array());

    BamRecordImplCoreTests::CheckRawData(bam1);
    BamRecordImplCoreTests::CheckRawData(bam2);
}

TEST(BAM_BamRecordImplCore, can_be_move_assigned)
{
    BamRecordImpl bam;
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    bam = std::move(BamRecordImplCoreTests::CreateBamImpl());
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    EXPECT_EQ(42, bam.Bin());
    EXPECT_EQ(42, bam.Flag());
    EXPECT_EQ(42, bam.InsertSize());
    EXPECT_EQ(42, bam.MapQuality());
    EXPECT_EQ(42, bam.MateReferenceId());
    EXPECT_EQ(42, bam.MatePosition());
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.ReferenceId());

    const TagCollection fetchedTags1 = bam.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    BamRecordImplCoreTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplCore, can_be_move_constructed)
{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    BamRecordImpl bam(std::move(BamRecordImplCoreTests::CreateBamImpl()));
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    EXPECT_EQ(42, bam.Bin());
    EXPECT_EQ(42, bam.Flag());
    EXPECT_EQ(42, bam.InsertSize());
    EXPECT_EQ(42, bam.MapQuality());
    EXPECT_EQ(42, bam.MateReferenceId());
    EXPECT_EQ(42, bam.MatePosition());
    EXPECT_EQ(42, bam.Position());
    EXPECT_EQ(42, bam.ReferenceId());

    const TagCollection fetchedTags1 = bam.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    BamRecordImplCoreTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplCore, can_set_and_query_alignment_flags)
{
    // same set of flags, different ways of getting there

    // raw number
    BamRecordImpl bam1;
    bam1.Flag(1107);

    // enum values
    BamRecordImpl bam2;
    bam2.Flag(BamRecordImpl::DUPLICATE | BamRecordImpl::MATE_1 | BamRecordImpl::REVERSE_STRAND |
              BamRecordImpl::PROPER_PAIR | BamRecordImpl::PAIRED);

    // convenience calls
    BamRecordImpl bam3;
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
