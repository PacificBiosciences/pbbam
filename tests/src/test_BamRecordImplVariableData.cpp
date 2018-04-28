// Author: Derek Barnett

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamRecordImpl.h>
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

    const auto rawData = PacBio::BAM::internal::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));

    EXPECT_EQ(expectedNameNulls, rawData->core.l_extranul);
    EXPECT_EQ(expectedNameLength, rawData->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps, rawData->core.n_cigar);
    EXPECT_EQ(expectedSeqLength, rawData->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, rawData->l_data);
}

}  // namespace BamRecordImplVariableDataTests

TEST(BamRecordImplVariableDataTest, InitEmpty)
{
    BamRecordImpl bam;
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, TagOnly_InitEmpty)
{
    BamRecordImpl bam;
    bam.Tags(TagCollection());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, TagOnly_InitNormal)
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

TEST(BamRecordImplVariableDataTest, TagOnly_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, TagOnly_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, TagOnly_ThenOverwriteWithEmptyTags)
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

TEST(BamRecordImplVariableDataTest, CigarOnly_InitEmpty)
{
    BamRecordImpl bam;
    bam.CigarData(std::string());
    EXPECT_EQ(0, bam.CigarData().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, CigarOnly_InitNormal_CigarObject)
{
    Cigar cigar;
    cigar.push_back(CigarOperation('=', 100));

    BamRecordImpl bam;
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData());
    EXPECT_TRUE("100=" == bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, CigarOnly_InitNormal_StdString)
{
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, CigarOnly_ThenOverwriteWithLongerCigar)
{
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, CigarOnly_ThenOverwriteWithShorterCigar)
{
    const std::string cigar = "100=";
    const std::string longerCigar = "100=10D100=10I100X";

    BamRecordImpl bam;
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, CigarOnly_ThenOverwriteWithEmptyCigar)
{
    const std::string cigar = "100=";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, CigarTag_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, CigarTag_Init_EmptyCigar)
{
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, CigarTag_Init_EmptyTag)
{
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, CigarTag_ThenOverwriteWithLongerCigar)
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

TEST(BamRecordImplVariableDataTest, CigarTag_ThenOverwriteWithShorterCigar)
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

TEST(BamRecordImplVariableDataTest, CigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, CigarTag_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, CigarTag_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, CigarTag_ThenOverwriteWithEmptyTags)
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

TEST(BamRecordImplVariableDataTest, SeqQualOnly_Init_Empty)
{
    BamRecordImpl bam;
    bam.SetSequenceAndQualities(std::string(), std::string());
    EXPECT_EQ(0, bam.Sequence().size());
    EXPECT_EQ(0, bam.Qualities().Fastq().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualOnly_Init_NormalQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualOnly_Init_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualOnly_Init_Preencoded)
{

    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    const size_t encodedLength = (sequence.size() + 1) / 2;
    char* encoded = static_cast<char*>(std::calloc(encodedLength, sizeof(char)));
    char* e = encoded;

    uint8_t nucleotideCode{};
    bool useHighWord = true;
    for (size_t i = 0; i < sequence.size(); ++i) {
        switch (sequence.at(i)) {
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

    if (encoded) free(encoded);
}

TEST(BamRecordImplVariableDataTest, SeqQualOnly_Init_Preencoded_EmptyQual)
{

    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";

    const auto encodedLength = (sequence.size() + 1) / 2;
    auto* encoded = static_cast<char*>(std::calloc(encodedLength, sizeof(char)));
    auto* e = encoded;

    uint8_t nucleotideCode{};
    bool useHighWord = true;
    for (size_t i = 0; i < sequence.size(); ++i) {
        switch (sequence.at(i)) {
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

    if (encoded) free(encoded);
}

TEST(BamRecordImplVariableDataTest, SeqQualOnly_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualOnly_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualOnly_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualOnly_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualOnly_ThenOverwriteWithEmptySeq)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualTag_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, SeqQualTag_Init_EmptySeqQual)
{
    const std::string sequence = "";
    const std::string qualities = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualTag_Init_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualTag_Init_EmptyTag)
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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithEmptySeq)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, SeqQualTag_ThenOverwriteWithEmptyTags)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigar_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigar_Init_EmptySeqQual)
{
    const std::string sequence = "";
    const std::string qualities = "";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualCigar_Init_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualCigar_Init_EmptyCigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithEmptySeq)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithLongerCigar)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithShorterCigar)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigar_ThenOverwriteWithEmptyCigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_Init_EmptySeqQual)
{
    const std::string sequence = "";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_Init_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_Init_EmptyCigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_Init_EmptyTag)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithEmptySeq)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerCigar)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterCigar)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, SeqQualCigarTag_ThenOverwriteWithEmptyTags)
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

TEST(BamRecordImplVariableDataTest, NameOnly_InitEmpty)
{
    BamRecordImpl bam;
    bam.Name(std::string());
    EXPECT_EQ(0, bam.Name().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameOnly_InitNormal)
{
    const std::string readName = "foo";

    BamRecordImpl bam;
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameOnly_ThenOverwriteWithLongerName)
{
    const std::string readName = "foo";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameOnly_ThenOverwriteWithShorterName)
{
    const std::string readName = "foo";
    const std::string longerName = "this is a long read name";

    BamRecordImpl bam;
    bam.Name(longerName);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameOnly_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string emptyName = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Name(emptyName);

    EXPECT_EQ(emptyName, bam.Name());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameTag_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, NameTag_Init_EmptyName)
{
    const std::string readName = "";

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

TEST(BamRecordImplVariableDataTest, NameTag_Init_EmptyTag)
{
    const std::string readName = "foo";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(0, bam.Tags().size());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameTag_ThenOverwriteWithLongerName)
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

TEST(BamRecordImplVariableDataTest, NameTag_ThenOverwriteWithShorterName)
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

TEST(BamRecordImplVariableDataTest, NameTag_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameTag_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, NameTag_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, NameTag_ThenOverwriteWithEmptyTags)
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

TEST(BamRecordImplVariableDataTest, NameCigar_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, NameCigar_Init_EmptyName)
{
    const std::string readName = "";
    const std::string cigar = "100=";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameCigar_Init_EmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameCigar_ThenOverwriteWithLongerName)
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

TEST(BamRecordImplVariableDataTest, NameCigar_ThenOverwriteWithShorterName)
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

TEST(BamRecordImplVariableDataTest, NameCigar_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameCigar_ThenOverwriteWithLongerCigar)
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

TEST(BamRecordImplVariableDataTest, NameCigar_ThenOverwriteWithShorterCigar)
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

TEST(BamRecordImplVariableDataTest, NameCigar_ThenOverwriteWithEmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameCigarTag_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_Init_EmptyName)
{
    const std::string readName = "";
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_Init_EmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar = "";

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

TEST(BamRecordImplVariableDataTest, NameCigarTag_Init_EmptyTag)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithLongerName)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithShorterName)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithLongerCigar)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithShorterCigar)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, NameCigarTag_ThenOverwriteWithEmptyTags)
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

TEST(BamRecordImplVariableDataTest, NameSeqQual_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, NameSeqQual_Init_EmptySeqQual)
{
    const std::string readName = "foo";
    const std::string sequence = "";
    const std::string qualities = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameSeqQual_Init_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithLongerName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithShorterName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(sequence, bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(shortSeq, bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameSeqQual_ThenOverwriteWithEmptySeq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty = "";

    BamRecordImpl bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities().Fastq());
    BamRecordImplVariableDataTests::CheckRawData(bam);
}

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_Init_EmptyName)
{
    const std::string readName = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_Init_EmptySeqQual)
{
    const std::string readName = "foo";
    const std::string sequence = "";
    const std::string qualities = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_Init_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_Init_EmptyTag)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithEmptySeq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualTag_ThenOverwriteWithEmptyTags)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_Init_EmptyName)
{
    const std::string readName = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_Init_EmptySeqQual)
{
    const std::string readName = "foo";
    const std::string sequence = "";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_Init_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_Init_EmptyCigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithEmptySeq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerCigar)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterCigar)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigar_ThenOverwriteWithEmptyCigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_Init_Normal)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_Init_EmptyName)
{
    const std::string readName = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_Init_EmptySeqQual)
{
    const std::string readName = "foo";
    const std::string sequence = "";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_Init_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_Init_EmptyCigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_Init_EmptyTag)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterName)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "";
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterSeq_NormalQual)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string shortSeq = "ACGT";
    const std::string shortQual = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptySeq)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerCigar)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterCigar)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string readName = "foo";
    const std::string sequence = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar = "100=";
    const std::string empty = "";

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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerTags)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterTags)
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

TEST(BamRecordImplVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptyTags)
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
