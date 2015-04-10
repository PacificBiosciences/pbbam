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

#ifdef PBBAM_TESTING
#define private public
#endif

#include <gtest/gtest.h>
#include <htslib/sam.h>
#include <pbbam/BamHeader.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace tests {

struct BamHdrDeleter
{
    void operator()(bam_hdr_t* hdr) {
        if (hdr)
            bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

} // namespace tests

TEST(BamHeaderTest, DefaultConstruction)
{
    BamHeader header;
    EXPECT_TRUE(header.Version().empty());
    EXPECT_TRUE(header.SortOrder().empty()); // default to unknown ?
    EXPECT_TRUE(header.ReadGroups().empty());
    EXPECT_TRUE(header.Sequences().empty());
    EXPECT_TRUE(header.Programs().empty());
    EXPECT_TRUE(header.Comments().empty());
}

TEST(BamHeaderTest, DecodeTest)
{
    const string& text = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                         "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                         "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                         "@RG\tID:rg1\tSM:control\n"
                         "@RG\tID:rg2\tSM:condition1\n"
                         "@RG\tID:rg3\tSM:condition1\n"
                         "@PG\tID:_foo_\tPN:ide\n"
                         "@CO\tipsum and so on\n"
                         "@CO\tcitation needed\n";

    BamHeader::SharedPtr header = BamHeader::FromSam(text);

    EXPECT_EQ(string("1.1"),       header->Version());
    EXPECT_EQ(string("queryname"), header->SortOrder());
    EXPECT_EQ(string("3.0b3"),     header->PacBioBamVersion());

    EXPECT_EQ(3, header->ReadGroups().size());
    EXPECT_TRUE(header->HasReadGroup("rg1"));
    EXPECT_TRUE(header->HasReadGroup("rg2"));
    EXPECT_TRUE(header->HasReadGroup("rg3"));

    EXPECT_EQ(string("control"),    header->ReadGroup("rg1").Sample());
    EXPECT_EQ(string("condition1"), header->ReadGroup("rg2").Sample());
    EXPECT_EQ(string("condition1"), header->ReadGroup("rg3").Sample());

    EXPECT_EQ(2, header->Sequences().size());
    EXPECT_TRUE(header->HasSequence("chr1"));
    EXPECT_TRUE(header->HasSequence("chr2"));
    EXPECT_EQ(string("chocobo"), header->Sequence("chr1").Species());
    EXPECT_EQ(string("chocobo"), header->Sequence("chr2").Species());
    EXPECT_EQ(string("2038"), header->Sequence("chr1").Length());
    EXPECT_EQ(string("3042"), header->Sequence("chr2").Length());

    EXPECT_EQ(1, header->Programs().size());
    EXPECT_TRUE(header->HasProgram("_foo_"));
    EXPECT_EQ(string("ide"), header->Program("_foo_").Name());

    EXPECT_EQ(2, header->Comments().size());
    EXPECT_EQ(string("ipsum and so on"), header->Comments().at(0));
    EXPECT_EQ(string("citation needed"), header->Comments().at(1));
}

TEST(BamHeaderCodecTest, EncodeTest)
{
    ReadGroupInfo rg1("rg1");
    rg1.Sample("control");
    ReadGroupInfo rg2("rg2");
    rg2.Sample("condition1");
    ReadGroupInfo rg3("rg3");
    rg3.Sample("condition1");

    SequenceInfo seq1("chr1");
    seq1.Length("2038").Species("chocobo");
    SequenceInfo seq2("chr2");
    seq2.Length("3042").Species("chocobo");

    ProgramInfo prog1("_foo_");
    prog1.Name("ide");

    BamHeader header;
    header.Version("1.1")
          .SortOrder("queryname")
          .PacBioBamVersion("3.0b3")
          .AddReadGroup(rg1)
          .AddReadGroup(rg2)
          .AddReadGroup(rg3)
          .AddSequence(seq1)
          .AddSequence(seq2)
          .AddProgram(prog1)
          .AddComment("ipsum and so on")
          .AddComment("citation needed");

    const string& expectedText = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                                 "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                                 "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                                 "@RG\tID:rg1\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:control\n"
                                 "@RG\tID:rg2\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                                 "@RG\tID:rg3\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                                 "@PG\tID:_foo_\tPN:ide\n"
                                 "@CO\tipsum and so on\n"
                                 "@CO\tcitation needed\n";

    const string& text = header.ToSam();
    EXPECT_EQ(expectedText, text);
}

TEST(BamHeaderTest, ConvertToRawDataOk)
{
    ReadGroupInfo rg1("rg1");
    rg1.Sample("control");
    ReadGroupInfo rg2("rg2");
    rg2.Sample("condition1");
    ReadGroupInfo rg3("rg3");
    rg3.Sample("condition1");

    SequenceInfo seq1("chr1");
    seq1.Length("2038").Species("chocobo");
    SequenceInfo seq2("chr2");
    seq2.Length("3042").Species("chocobo");

    ProgramInfo prog1("_foo_");
    prog1.Name("ide");

    BamHeader header;
    header.Version("1.1")
          .SortOrder("queryname")
          .PacBioBamVersion("3.0b3")
          .AddReadGroup(rg1)
          .AddReadGroup(rg2)
          .AddReadGroup(rg3)
          .AddSequence(seq1)
          .AddSequence(seq2)
          .AddProgram(prog1)
          .AddComment("ipsum and so on")
          .AddComment("citation needed");

    const string& expectedText = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                                 "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                                 "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                                 "@RG\tID:rg1\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:control\n"
                                 "@RG\tID:rg2\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                                 "@RG\tID:rg3\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                                 "@PG\tID:_foo_\tPN:ide\n"
                                 "@CO\tipsum and so on\n"
                                 "@CO\tcitation needed\n";


    const string& text = header.ToSam();
    PBBAM_SHARED_PTR<bam_hdr_t> rawData(sam_hdr_parse(text.size(), text.c_str()), tests::BamHdrDeleter());
    rawData->ignore_sam_err = 0;
    rawData->cigar_tab = NULL;
    rawData->l_text = text.size();
    rawData->text = (char*)calloc(rawData->l_text + 1, 1);
    memcpy(rawData->text, text.c_str(), rawData->l_text);

    const string& rawText = string(rawData->text, rawData->l_text);
    EXPECT_EQ(expectedText, rawText);
}

TEST(BamHeaderTest, ExtractFromRawDataOk)
{
    ReadGroupInfo rg1("rg1");
    rg1.Sample("control");
    ReadGroupInfo rg2("rg2");
    rg2.Sample("condition1");
    ReadGroupInfo rg3("rg3");
    rg3.Sample("condition1");

    SequenceInfo seq1("chr1");
    seq1.Length("2038").Species("chocobo");
    SequenceInfo seq2("chr2");
    seq2.Length("3042").Species("chocobo");

    ProgramInfo prog1("_foo_");
    prog1.Name("ide");

    BamHeader header;
    header.Version("1.1")
          .SortOrder("queryname")
          .PacBioBamVersion("3.0b3")
          .AddReadGroup(rg1)
          .AddReadGroup(rg2)
          .AddReadGroup(rg3)
          .AddSequence(seq1)
          .AddSequence(seq2)
          .AddProgram(prog1)
          .AddComment("ipsum and so on")
          .AddComment("citation needed");

    const string& expectedText = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                                 "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                                 "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                                 "@RG\tID:rg1\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:control\n"
                                 "@RG\tID:rg2\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                                 "@RG\tID:rg3\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                                 "@PG\tID:_foo_\tPN:ide\n"
                                 "@CO\tipsum and so on\n"
                                 "@CO\tcitation needed\n";


    string text = header.ToSam();
    PBBAM_SHARED_PTR<bam_hdr_t> rawData(sam_hdr_parse(text.size(), text.c_str()), tests::BamHdrDeleter());
    rawData->ignore_sam_err = 0;
    rawData->cigar_tab = NULL;
    rawData->l_text = text.size();
    rawData->text = (char*)calloc(rawData->l_text + 1, 1);
    memcpy(rawData->text, text.c_str(), rawData->l_text);

    const BamHeader::SharedPtr& newHeader = BamHeader::FromSam(string(rawData->text, rawData->l_text));

    EXPECT_EQ(header.Version(),          newHeader->Version());
    EXPECT_EQ(header.SortOrder(),        newHeader->SortOrder());
    EXPECT_EQ(header.PacBioBamVersion(), newHeader->PacBioBamVersion());

    text = newHeader->ToSam();
    EXPECT_EQ(expectedText, text);
}
