// Author: Derek Barnett

#include <pbbam/BamRecordImpl.h>

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamTagCodec.h>
#include <pbbam/SamTagCodec.h>
#include <pbbam/Tag.h>
#include <pbbam/TagCollection.h>

#include "../src/MemoryUtils.h"

using namespace PacBio;
using namespace PacBio::BAM;

// NOTE: this file has a *TON* of tests. Probably overkill, but I wanted to check
//       every possible combination of variable data, and then manipulate each
//       element within each combo to shrink & expand.

namespace BamRecordImplVariableDataTests {

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
    // l_qname + (n_cigar * 4) + (l_qseq+1)/2 + l_qseq + <encoded length>
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

}  // namespace BamRecordImplVariableDataTests

TEST(BAM_BamRecordImplVarData, default_construction_all_data_is_empty_and_valid)
{
    BamRecordImpl bam;
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, default_construction_tags_are_empty)
{
    BamRecordImpl bam;
    bam.Tags(TagCollection());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_add_tags)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(tags);

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);
}

TEST(BAM_BamRecordImplVarData, can_add_tags_then_overwrite_with_longer_tags)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(tags);
    bam.Tags(longerTags);

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);
}

TEST(BAM_BamRecordImplVarData, can_add_tags_then_overwrite_with_shorter_tags)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(longerTags);
    bam.Tags(tags);

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);
}

TEST(BAM_BamRecordImplVarData, can_add_tags_then_overwrite_with_empty_tags)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, default_construction_cigar_is_empty)
{
    BamRecordImpl bam;
    bam.CigarData(std::string());
    EXPECT_EQ(0, bam.CigarData().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_from_cigar_object)
{
    Cigar cigar;
    cigar.push_back(CigarOperation('=', 100));

    BamRecordImpl bam;
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData());
    EXPECT_TRUE("100=" == bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_from_cigar_string)
{
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_then_overwrite_with_longer_cigar)
{
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_then_overwrite_with_shorter_cigar)
{
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_then_overwrite_with_empty_cigar)
{
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag)
{
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_with_empty_cigar_string)
{
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(empty, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_with_empty_cigar)
{
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_then_overwrite_with_longer_cigar)
{
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(longerCigar);

    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_then_overwrite_with_shorter_cigar)
{
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.CigarData(longerCigar);
    bam.Tags(tags);
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_then_overwrite_with_empty_cigar)
{
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(empty, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_then_overwrite_with_longer_tag_collection)
{
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_then_overwrite_with_shorter_tag_collection)
{
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_cigar_tag_then_overwrite_with_empty_tag_collection)
{
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_with_both_empty)
{
    BamRecordImpl bam;
    bam.SetSequenceAndQualities(std::string(), std::string());
    EXPECT_EQ(0, bam.Sequence().size());
    EXPECT_EQ(0, bam.Qualities().Fastq().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_with_both_normal)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_with_empty_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_preencoded)
{

    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    const size_t encodedLength = (sequence.size() + 1) / 2;
    char* encoded = static_cast<char*>(std::calloc(encodedLength, sizeof(char)));
    char* e = encoded;

    uint8_t nucleotideCode{};
    bool useHighWord = true;
    for (char i : sequence) {
        switch (i) {
            case 'A':
                nucleotideCode = 1;
                break;
            case 'C':
                nucleotideCode = 2;
                break;
            case 'G':
                nucleotideCode = 4;
                break;
            case 'T':
                nucleotideCode = 8;
                break;
            default:
                EXPECT_FALSE(true);
                break;
        }

        // pack the nucleotide code
        if (useHighWord) {
            *e = nucleotideCode << 4;
            useHighWord = false;
        } else {
            *e |= nucleotideCode;
            ++e;
            useHighWord = true;
        }
    }

    BamRecordImpl bam;
    bam.SetPreencodedSequenceAndQualities(encoded, sequence.size(), qualities.c_str());

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);

    free(encoded);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_preencoded_with_empty_qual)
{

    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;

    const auto encodedLength = (sequence.size() + 1) / 2;
    auto* encoded = static_cast<char*>(std::calloc(encodedLength, sizeof(char)));
    auto* e = encoded;

    uint8_t nucleotideCode{};
    bool useHighWord = true;
    for (char i : sequence) {
        switch (i) {
            case 'A':
                nucleotideCode = 1;
                break;
            case 'C':
                nucleotideCode = 2;
                break;
            case 'G':
                nucleotideCode = 4;
                break;
            case 'T':
                nucleotideCode = 8;
                break;
            default:
                EXPECT_FALSE(true);
                break;
        }

        // pack the nucleotide code
        if (useHighWord) {
            *e = nucleotideCode << 4;
            useHighWord = false;
        } else {
            *e |= nucleotideCode;
            ++e;
            useHighWord = true;
        }
    }

    BamRecordImpl bam;
    bam.SetPreencodedSequenceAndQualities(encoded, sequence.size(), qualities.c_str());

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);

    free(encoded);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_overwrite_with_longer_seq_and_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_overwrite_with_longer_seq_only)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_overwrite_with_shorter_seq_and_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_overwrite_with_shorter_seq_only)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_overwrite_with_empty_seq_)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty;

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_seq_qual_then_tags)
{
    const std::string sequence;
    const std::string qualities;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_with_empty_qual_then_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_empty_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_longer_seq_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_longer_seq_qual_with_empty_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_shorter_seq_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_shorter_seq_qual_with_empty_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_empty_seq)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_longer_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_shorter_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_tags_then_empty_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_seq_qual_then_cigar)
{
    const std::string sequence;
    const std::string qualities;
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_and_empty_qual_then_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_empty_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar;

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_overwrite_with_longer_seq_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_overwrite_with_longer_seq_empty_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_overwrite_with_shorter_seq_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_overwrite_with_shorter_seq_empty_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_overwrite_with_empty_seq)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_overwrite_with_longer_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_overwrite_with_shorter_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_overwrite_with_empty_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_with_empty_seq_qual)
{
    const std::string sequence;
    const std::string qualities;
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_with_empty_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_with_empty_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_with_empty_tag)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_seq_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_seq_only)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_seq_qual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_seq_only)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_empty_seq)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(longerCigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.Tags(tags);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_empty_cigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_seq_qual_then_cigar_then_tag_then_overwrite_with_empty_tags)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_with_empty_name)
{
    BamRecordImpl bam;
    bam.Name(std::string());
    EXPECT_EQ(0, bam.Name().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_with_name_only)
{
    const std::string readName = "foo";

    BamRecordImpl bam;
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string emptyName;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Name(emptyName);

    EXPECT_EQ(emptyName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_tag)
{
    const std::string readName = "foo";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_name_then_tag)
{
    const std::string readName;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_tag)
{
    const std::string readName = "foo";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_tag_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(tags);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_tag_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.Tags(tags);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_tag_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(tags);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_tag_then_overwrite_with_longer_tags)
{
    const std::string readName = "foo";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(readName, bam.Name());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_tag_then_overwrite_with_shorter_tags)
{
    const std::string readName = "foo";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_tag_then_overwrite_with_empty_tags)
{
    const std::string readName = "foo";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_name_then_cigar)
{
    const std::string readName;
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_name_then_empty_cigar)
{
    const std::string readName = "foo";
    const std::string cigar;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.CigarData(cigar);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_overwrite_with_longer_cigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_overwrite_with_shorter_cigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_overwrite_with_empty_cigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_name_then_cigar_then_tag)
{
    const std::string readName;
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_cigar_then_tag)
{
    const std::string readName = "foo";
    const std::string cigar;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_empty_tag)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_longer_cigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_shorter_cigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(longerCigar);
    bam.Tags(tags);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_empty_cigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_longer_tag)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_shorter_tag)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_cigar_then_tag_then_overwrite_with_empty_tag)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence;
    const std::string qualities;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_longer_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_longer_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_shorter_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_shorter_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_overwrite_with_empty_seq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_name_then_seq_qual_then_tag)
{
    const std::string readName;
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_seq_qual_then_tag)
{
    const std::string readName = "foo";
    const std::string sequence;
    const std::string qualities;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_only_then_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_empty_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(0, bam.Tags().size());

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_tag_then_overwrite_with_longer_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_tag_then_overwrite_with_longer_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_tag_then_overwrite_with_shorter_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_tag_then_overwrite_with_shorter_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag_then_overwrite_with_empty_seq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag_then_overwrite_with_longer_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag_then_overwrite_with_shorter_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_tag_then_overwrite_with_empty_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_name_then_seq_qual_then_cigar)
{
    const std::string readName;
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_seq_qual_then_cigar)
{
    const std::string readName = "foo";
    const std::string sequence;
    const std::string qualities;
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_only_then_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_empty_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_cigar_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_longer_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_longer_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_shorter_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_shorter_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_cigar_then_overwrite_with_empty_seq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_longer_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_shorter_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_overwrite_with_empty_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_cigar_then_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_empty_name_then_seq_qual_then_cigar_then_tag)
{
    const std::string readName;
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_seq_qual_then_cigar_then_tag)
{
    const std::string readName = "foo";
    const std::string sequence;
    const std::string qualities;
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_empty_seq_only_then_cigar_then_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_empty_cigar_then_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData, can_set_name_then_seq_qual_then_cigar_then_empty_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_empty_name)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities;
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_seq_qual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_seq_only)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_empty_seq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.Tags(tags);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_empty_cigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_longer_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_shorter_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = int32_t{-42};

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BAM_BamRecordImplVarData,
     can_set_name_then_seq_qual_then_cigar_then_tag_then_overwrite_with_empty_tag)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}
