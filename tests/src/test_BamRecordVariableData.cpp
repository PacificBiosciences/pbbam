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

#include <gtest/gtest.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamTagCodec.h>
#include <pbbam/SamTagCodec.h>
#include <pbbam/Tag.h>
#include <pbbam/TagCollection.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;

// NOTE: this file has a *TON* of tests. Probably overkill, but I wanted to check
//       every possible combination of variable data, and then manipulate each
//       element within each combo to shrink & expand.

namespace tests {

static
void CheckRawData(const BamRecord& bam)
{
    // ensure raw data (lengths at least) matches API-facing data

    const uint32_t expectedNameLength  = bam.Name().size() + 1;
    const uint32_t expectedNumCigarOps = bam.CigarData().size();
    const int32_t  expectedSeqLength   = bam.Sequence().length();
    const size_t   expectedTagsLength  = BamTagCodec::Encode(bam.Tags()).size();

    //  Name        CIGAR         Sequence       Quals      Tags
    // l_qname + (n_cigar * 4) + (l_qseq+1)/2 + l_qseq + <encoded length>

    const int expectedTotalDataLength = expectedNameLength +
                                        (expectedNumCigarOps * 4) +
                                        (expectedSeqLength+1)/2 +
                                         expectedSeqLength +
                                         expectedTagsLength;

    EXPECT_EQ(expectedNameLength,      bam.RawData()->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps,     bam.RawData()->core.n_cigar);
    EXPECT_EQ(expectedSeqLength,       bam.RawData()->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, bam.RawData()->l_data);
}

} // namespace tests

TEST(BamRecordVariableDataTest, InitEmpty)
{
    BamRecord bam;
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, TagOnly_InitEmpty)
{
    BamRecord bam;
    bam.Tags(TagCollection());
    EXPECT_EQ(0, bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, TagOnly_InitNormal)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

TEST(BamRecordVariableDataTest, TagOnly_ThenOverwriteWithLongerTags)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
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

TEST(BamRecordVariableDataTest, TagOnly_ThenOverwriteWithShorterTags)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Tags(longerTags);
    bam.Tags(tags);

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);
}

TEST(BamRecordVariableDataTest, TagOnly_ThenOverwriteWithEmptyTags)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(0, bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarOnly_InitEmpty)
{
    BamRecord bam;
    bam.CigarData(std::string());
    EXPECT_EQ(0, bam.CigarData().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarOnly_InitNormal_CigarObject)
{
    Cigar cigar;
    cigar.push_back(CigarOperation('M', 100));

    BamRecord bam;
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData());
    EXPECT_TRUE("100M" == bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarOnly_InitNormal_StdString)
{
    const std::string cigar = "100M";

    BamRecord bam;
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarOnly_ThenOverwriteWithLongerCigar)
{
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarOnly_ThenOverwriteWithShorterCigar)
{
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarOnly_ThenOverwriteWithEmptyCigar)
{
    const std::string cigar = "100M";
    const std::string empty = "";

    BamRecord bam;
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(empty, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_Init_Normal)
{
    const std::string cigar = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_Init_EmptyCigar)
{
    const std::string cigar = "100M";
    const std::string empty = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_Init_EmptyTag)
{
    const std::string cigar = "100M";

    BamRecord bam;
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0,     bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_ThenOverwriteWithLongerCigar)
{
    const std::string cigar = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_ThenOverwriteWithShorterCigar)
{
    const std::string cigar = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string cigar = "100M";
    const std::string empty = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_ThenOverwriteWithLongerTags)
{
    const std::string cigar = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_ThenOverwriteWithShorterTags)
{
    const std::string cigar = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, CigarTag_ThenOverwriteWithEmptyTags)
{
    const std::string cigar = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    EXPECT_EQ(0,     bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_Init_Empty)
{
    BamRecord bam;
    bam.SetSequenceAndQualities(std::string(), std::string());
    EXPECT_EQ(0, bam.Sequence().size());
    EXPECT_EQ(0, bam.Qualities().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_Init_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_Init_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_Init_Preencoded) {

    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    const size_t encodedLength = static_cast<size_t>((sequence.size()+1)/2);
    char* encoded = (char*)::calloc(encodedLength, sizeof(char));
    char* e = encoded;

    uint8_t nucleotideCode;
    bool useHighWord = true;
    for (size_t i = 0; i < sequence.size(); ++i) {
        switch (sequence.at(i)) {
            case 'A' : nucleotideCode = 1;  break;
            case 'C' : nucleotideCode = 2;  break;
            case 'G' : nucleotideCode = 4;  break;
            case 'T' : nucleotideCode = 8;  break;
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

    BamRecord bam;
    bam.SetPreencodedSequenceAndQualities(encoded, sequence.size(), qualities.c_str());

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_Init_Preencoded_EmptyQual) {

    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";

    const size_t encodedLength = static_cast<size_t>((sequence.size()+1)/2);
    char* encoded = (char*)::calloc(encodedLength, sizeof(char));
    char* e = encoded;

    uint8_t nucleotideCode;
    bool useHighWord = true;
    for (size_t i = 0; i < sequence.size(); ++i) {
        switch (sequence.at(i)) {
            case 'A' : nucleotideCode = 1;  break;
            case 'C' : nucleotideCode = 2;  break;
            case 'G' : nucleotideCode = 4;  break;
            case 'T' : nucleotideCode = 8;  break;
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

    BamRecord bam;
    bam.SetPreencodedSequenceAndQualities(encoded, sequence.size(), qualities.c_str());

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualOnly_ThenOverwriteWithEmptySeq)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty     = "";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_Init_Normal)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_Init_EmptySeqQual)
{
    const std::string sequence  = "";
    const std::string qualities = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_Init_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_Init_EmptyTag)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(0,         bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithEmptySeq)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithLongerTags)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithShorterTags)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualTag_ThenOverwriteWithEmptyTags)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(0,         bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_Init_Normal)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_Init_EmptySeqQual)
{
    const std::string sequence  = "";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_Init_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_Init_EmptyCigar)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithEmptySeq)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithLongerCigar)
{
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(sequence,    bam.Sequence());
    EXPECT_EQ(qualities,   bam.Qualities());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithShorterCigar)
{
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigar_ThenOverwriteWithEmptyCigar)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(empty,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_Init_Normal)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_Init_EmptySeqQual)
{
    const std::string sequence  = "";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_Init_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_Init_EmptyCigar)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_Init_EmptyTag)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(0,         bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithEmptySeq)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(empty, bam.Sequence());
    EXPECT_EQ(empty, bam.Qualities());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerCigar)
{
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(longerCigar);

    EXPECT_EQ(sequence,    bam.Sequence());
    EXPECT_EQ(qualities,   bam.Qualities());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterCigar)
{
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.Tags(tags);
    bam.CigarData(cigar);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(empty,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithLongerTags)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithShorterTags)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, SeqQualCigarTag_ThenOverwriteWithEmptyTags)
{
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(0,         bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameOnly_InitEmpty)
{
    BamRecord bam;
    bam.Name(std::string());
    EXPECT_EQ(0, bam.Name().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameOnly_InitNormal)
{
    const std::string readName = "foo";

    BamRecord bam;
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameOnly_ThenOverwriteWithLongerName)
{
    const std::string readName   = "foo";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(readName);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameOnly_ThenOverwriteWithShorterName)
{
    const std::string readName   = "foo";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(longerName);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameOnly_ThenOverwriteWithEmptyName)
{
    const std::string readName  = "foo";
    const std::string emptyName = "";

    BamRecord bam;
    bam.Name(readName);
    bam.Name(emptyName);

    EXPECT_EQ(emptyName, bam.Name());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_Init_Normal)
{
    const std::string readName = "foo";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_Init_EmptyName)
{
    const std::string readName = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_Init_EmptyTag)
{
    const std::string readName = "foo";

    BamRecord bam;
    bam.Name(readName);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(0,        bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_ThenOverwriteWithLongerName)
{
    const std::string readName   = "foo";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_ThenOverwriteWithShorterName)
{
    const std::string readName   = "foo";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string empty    = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_ThenOverwriteWithLongerTags)
{
    const std::string readName = "foo";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_ThenOverwriteWithShorterTags)
{
    const std::string readName = "foo";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameTag_ThenOverwriteWithEmptyTags)
{
    const std::string readName = "foo";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(0,        bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_Init_Normal)
{
    const std::string readName = "foo";
    const std::string cigar   = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_Init_EmptyName)
{
    const std::string readName = "";
    const std::string cigar    = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_Init_EmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar    = "";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_ThenOverwriteWithLongerName)
{
    const std::string readName   = "foo";
    const std::string cigar      = "100M";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(cigar,      bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_ThenOverwriteWithShorterName)
{
    const std::string readName   = "foo";
    const std::string cigar      = "100M";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(longerName);
    bam.CigarData(cigar);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";
    const std::string empty    = "";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Name(empty);

    EXPECT_EQ(empty, bam.Name());
    EXPECT_EQ(cigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_ThenOverwriteWithLongerCigar)
{
    const std::string readName    = "foo";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName,    bam.Name());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_ThenOverwriteWithShorterCigar)
{
    const std::string readName    = "foo";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigar_ThenOverwriteWithEmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";
    const std::string empty    = "";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty,    bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_Init_Normal)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_Init_EmptyName)
{
    const std::string readName = "";
    const std::string cigar    = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_Init_EmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar    = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_Init_EmptyTag)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    EXPECT_EQ(0,        bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithLongerName)
{
    const std::string readName   = "foo";
    const std::string cigar      = "100M";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(cigar,      bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithShorterName)
{
    const std::string readName   = "foo";
    const std::string cigar      = "100M";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(longerName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(readName);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithEmptyName)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";
    const std::string empty    = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
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

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithLongerCigar)
{
    const std::string readName    = "foo";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName,    bam.Name());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithShorterCigar)
{
    const std::string readName    = "foo";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(longerCigar);
    bam.Tags(tags);
    bam.CigarData(cigar);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";
    const std::string empty    = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithLongerTags)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithShorterTags)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameCigarTag_ThenOverwriteWithEmptyTags)
{
    const std::string readName = "foo";
    const std::string cigar    = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    EXPECT_EQ(0,        bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_Init_Normal)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_Init_EmptySeqQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "";
    const std::string qualities = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_Init_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithLongerName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(sequence,   bam.Sequence());
    EXPECT_EQ(qualities,  bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithShorterName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Name(readName);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithEmptyName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty     = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Name(empty);

    EXPECT_EQ(empty,     bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQual_ThenOverwriteWithEmptySeq)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty     = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty,    bam.Sequence());
    EXPECT_EQ(empty,    bam.Qualities());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_Init_Normal)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_Init_EmptyName)
{
    const std::string readName  = "";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_Init_EmptySeqQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "";
    const std::string qualities = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_Init_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_Init_EmptyTag)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(0,         bam.Tags().size());

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(sequence,   bam.Sequence());
    EXPECT_EQ(qualities,  bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Name(readName);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithEmptyName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Name(empty);

    EXPECT_EQ(empty,     bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithEmptySeq)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty,    bam.Sequence());
    EXPECT_EQ(empty,    bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithLongerTags)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithShorterTags)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualTag_ThenOverwriteWithEmptyTags)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(0,         bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_Init_Normal)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_Init_EmptyName)
{
    const std::string readName  = "";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_Init_EmptySeqQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_Init_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_Init_EmptyCigar)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerName)
{
    const std::string readName   = "foo";
    const std::string sequence   = "ACGTACGTACGT";
    const std::string qualities  = "?]?]?]?]?]?]";
    const std::string cigar      = "100M";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Name(longerName);

    EXPECT_EQ(longerName,  bam.Name());
    EXPECT_EQ(sequence,    bam.Sequence());
    EXPECT_EQ(qualities,   bam.Qualities());
    EXPECT_EQ(cigar,       bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterName)
{
    const std::string readName   = "foo";
    const std::string sequence   = "ACGTACGTACGT";
    const std::string qualities  = "?]?]?]?]?]?]";
    const std::string cigar      = "100M";
    const std::string longerName = "this is a long read name";

    BamRecord bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Name(readName);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,    bam.Sequence());
    EXPECT_EQ(qualities,   bam.Qualities());
    EXPECT_EQ(cigar,       bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithEmptyName)
{
    const std::string readName   = "foo";
    const std::string sequence   = "ACGTACGTACGT";
    const std::string qualities  = "?]?]?]?]?]?]";
    const std::string cigar      = "100M";
    const std::string empty      = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Name(empty);

    EXPECT_EQ(empty,       bam.Name());
    EXPECT_EQ(sequence,    bam.Sequence());
    EXPECT_EQ(qualities,   bam.Qualities());
    EXPECT_EQ(cigar,       bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string readName   = "foo";
    const std::string sequence   = "ACGTACGTACGT";
    const std::string qualities  = "?]?]?]?]?]?]";
    const std::string cigar      = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName   = "foo";
    const std::string sequence   = "ACGTACGTACGT";
    const std::string qualities  = "";
    const std::string cigar      = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string readName   = "foo";
    const std::string sequence   = "ACGTACGTACGT";
    const std::string qualities  = "?]?]?]?]?]?]";
    const std::string cigar      = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName   = "foo";
    const std::string sequence   = "ACGTACGTACGT";
    const std::string qualities  = "?]?]?]?]?]?]";
    const std::string cigar      = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithEmptySeq)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty,    bam.Sequence());
    EXPECT_EQ(empty,    bam.Qualities());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithLongerCigar)
{
    const std::string readName    = "foo";
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName,    bam.Name());
    EXPECT_EQ(sequence,    bam.Sequence());
    EXPECT_EQ(qualities,   bam.Qualities());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithShorterCigar)
{
    const std::string readName    = "foo";
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.CigarData(cigar);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigar_ThenOverwriteWithEmptyCigar)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.CigarData(empty);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(empty,     bam.CigarData().ToStdString());
    tests::CheckRawData(bam);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_Init_Normal)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_Init_EmptyName)
{
    const std::string readName  = "";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_Init_EmptySeqQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_Init_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_Init_EmptyCigar)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_Init_EmptyTag)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(0,         bam.Tags().size());
    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(longerName);

    EXPECT_EQ(longerName, bam.Name());
    EXPECT_EQ(sequence,   bam.Sequence());
    EXPECT_EQ(qualities,  bam.Qualities());
    EXPECT_EQ(cigar,      bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string longerName = "this is a long read name";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(longerName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(readName);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptyName)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Name(empty);

    EXPECT_EQ(empty,     bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerSeq_NormalQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerSeq_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(shortSeq, shortQual);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(sequence, qualities);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterSeq_NormalQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "?]?]";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterSeq_EmptyQual)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string shortSeq  = "ACGT";
    const std::string shortQual = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(shortSeq, shortQual);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(shortSeq,  bam.Sequence());
    EXPECT_EQ(shortQual, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptySeq)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.SetSequenceAndQualities(empty, empty);

    EXPECT_EQ(readName, bam.Name());
    EXPECT_EQ(empty,    bam.Sequence());
    EXPECT_EQ(empty,    bam.Qualities());
    EXPECT_EQ(cigar,    bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerCigar)
{
    const std::string readName    = "foo";
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(longerCigar);

    EXPECT_EQ(readName,    bam.Name());
    EXPECT_EQ(sequence,    bam.Sequence());
    EXPECT_EQ(qualities,   bam.Qualities());
    EXPECT_EQ(longerCigar, bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterCigar)
{
    const std::string readName    = "foo";
    const std::string sequence    = "ACGTACGTACGT";
    const std::string qualities   = "?]?]?]?]?]?]";
    const std::string cigar       = "100M";
    const std::string longerCigar = "100=10D100M10I100X";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(longerCigar);
    bam.Tags(tags);
    bam.CigarData(cigar);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptyCigar)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";
    const std::string empty     = "";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.CigarData(empty);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(empty,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithLongerTags)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(longerTags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithShorterTags)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection longerTags;
    longerTags["HX"] = std::string("1abc75");
    longerTags["HX"].Modifier(TagModifier::HEX_STRING);
    longerTags["CA"] = std::vector<uint8_t>({34, 5, 125});
    longerTags["XY"] = (int32_t)-42;

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(longerTags);
    bam.Tags(tags);

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());

    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");

    const std::string sam = SamTagCodec::Encode(bam.Tags());
    EXPECT_EQ(expected, sam);

    tests::CheckRawData(bam);
}

TEST(BamRecordVariableDataTest, NameSeqQualCigarTag_ThenOverwriteWithEmptyTags)
{
    const std::string readName  = "foo";
    const std::string sequence  = "ACGTACGTACGT";
    const std::string qualities = "?]?]?]?]?]?]";
    const std::string cigar     = "100M";

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecord bam;
    bam.Name(readName);
    bam.SetSequenceAndQualities(sequence, qualities);
    bam.CigarData(cigar);
    bam.Tags(tags);
    bam.Tags(TagCollection());

    EXPECT_EQ(readName,  bam.Name());
    EXPECT_EQ(sequence,  bam.Sequence());
    EXPECT_EQ(qualities, bam.Qualities());
    EXPECT_EQ(cigar,     bam.CigarData().ToStdString());
    EXPECT_EQ(0,         bam.Tags().size());
    tests::CheckRawData(bam);
}
