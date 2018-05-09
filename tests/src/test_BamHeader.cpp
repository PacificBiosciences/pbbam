// Author: Derek Barnett

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include <htslib/sam.h>

#include <pbbam/BamHeader.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace BamHeaderTests {

struct BamHdrDeleter
{
    void operator()(bam_hdr_t* hdr)
    {
        if (hdr) bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

}  // namespace BamHeaderTests

TEST(BamHeaderTest, DefaultConstruction)
{
    BamHeader header;
    EXPECT_TRUE(header.Version().empty());
    EXPECT_TRUE(header.SortOrder().empty());  // default to unknown ?
    EXPECT_TRUE(header.ReadGroups().empty());
    EXPECT_TRUE(header.Sequences().empty());
    EXPECT_TRUE(header.Programs().empty());
    EXPECT_TRUE(header.Comments().empty());

    EXPECT_THROW(header.Program("foo"), std::exception);
    EXPECT_THROW(header.ReadGroup("foo"), std::exception);
    EXPECT_THROW(header.SequenceId("foo"), std::exception);
    EXPECT_THROW(header.SequenceLength(42), std::exception);
    EXPECT_THROW(header.SequenceName(42), std::exception);
}

TEST(BamHeaderTest, DecodeTest)
{
    const std::string text{
        "@HD\tVN:1.1\tSO:queryname\tpb:3.0.1\n"
        "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
        "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
        "@RG\tID:rg1\tSM:control\n"
        "@RG\tID:rg2\tSM:condition1\n"
        "@RG\tID:rg3\tSM:condition1\n"
        "@PG\tID:_foo_\tPN:ide\n"
        "@CO\tipsum and so on\n"
        "@CO\tcitation needed\n"};

    BamHeader header = BamHeader(text);

    EXPECT_EQ(std::string("1.1"), header.Version());
    EXPECT_EQ(std::string("queryname"), header.SortOrder());
    EXPECT_EQ(std::string("3.0.1"), header.PacBioBamVersion());

    EXPECT_EQ(3, header.ReadGroups().size());
    EXPECT_TRUE(header.HasReadGroup("rg1"));
    EXPECT_TRUE(header.HasReadGroup("rg2"));
    EXPECT_TRUE(header.HasReadGroup("rg3"));

    EXPECT_EQ(std::string("control"), header.ReadGroup("rg1").Sample());
    EXPECT_EQ(std::string("condition1"), header.ReadGroup("rg2").Sample());
    EXPECT_EQ(std::string("condition1"), header.ReadGroup("rg3").Sample());

    EXPECT_EQ(2, header.Sequences().size());
    EXPECT_TRUE(header.HasSequence("chr1"));
    EXPECT_TRUE(header.HasSequence("chr2"));
    EXPECT_EQ(std::string("chocobo"), header.Sequence("chr1").Species());
    EXPECT_EQ(std::string("chocobo"), header.Sequence("chr2").Species());
    EXPECT_EQ(std::string("2038"), header.Sequence("chr1").Length());
    EXPECT_EQ(std::string("3042"), header.Sequence("chr2").Length());

    EXPECT_EQ(1, header.Programs().size());
    EXPECT_TRUE(header.HasProgram("_foo_"));
    EXPECT_EQ(std::string("ide"), header.Program("_foo_").Name());

    EXPECT_EQ(2, header.Comments().size());
    EXPECT_EQ(std::string("ipsum and so on"), header.Comments().at(0));
    EXPECT_EQ(std::string("citation needed"), header.Comments().at(1));
}

TEST(BamHeaderTest, VersionCheckOk)
{
    auto expectFail = [](std::string&& label, std::string&& text) {
        SCOPED_TRACE(label);
        EXPECT_THROW(BamHeader{text}, std::runtime_error);
    };
    expectFail("empty version", "@HD\tVN:1.1\tSO:queryname\tpb:\n");
    expectFail("old beta version", "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n");
    expectFail("old beta version", "@HD\tVN:1.1\tSO:queryname\tpb:3.0b7\n");
    expectFail("invalid value", "@HD\tVN:1.1\tSO:queryname\tpb:3.0.should_not_work\n");
    expectFail("earlier than minimum", "@HD\tVN:1.1\tSO:queryname\tpb:3.0.0\n");

    // correct version syntax, number
    EXPECT_NO_THROW(BamHeader{"@HD\tVN:1.1\tSO:queryname\tpb:3.0.1\n"});
}

TEST(BamHeaderTest, EncodeTest)
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
        .PacBioBamVersion("3.0.1")
        .AddReadGroup(rg1)
        .AddReadGroup(rg2)
        .AddReadGroup(rg3)
        .AddSequence(seq1)
        .AddSequence(seq2)
        .AddProgram(prog1)
        .AddComment("ipsum and so on")
        .AddComment("citation needed");

    const std::string expectedText{
        "@HD\tVN:1.1\tSO:queryname\tpb:3.0.1\n"
        "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
        "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
        "@RG\tID:rg1\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:control\tPM:SEQUEL\n"
        "@RG\tID:rg2\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\tPM:SEQUEL\n"
        "@RG\tID:rg3\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\tPM:SEQUEL\n"
        "@PG\tID:_foo_\tPN:ide\n"
        "@CO\tipsum and so on\n"
        "@CO\tcitation needed\n"};

    const std::string text = header.ToSam();
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
        .PacBioBamVersion("3.0.1")
        .AddReadGroup(rg1)
        .AddReadGroup(rg2)
        .AddReadGroup(rg3)
        .AddSequence(seq1)
        .AddSequence(seq2)
        .AddProgram(prog1)
        .AddComment("ipsum and so on")
        .AddComment("citation needed");

    const std::string expectedText{
        "@HD\tVN:1.1\tSO:queryname\tpb:3.0.1\n"
        "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
        "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
        "@RG\tID:rg1\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:control\tPM:SEQUEL\n"
        "@RG\tID:rg2\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\tPM:SEQUEL\n"
        "@RG\tID:rg3\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\tPM:SEQUEL\n"
        "@PG\tID:_foo_\tPN:ide\n"
        "@CO\tipsum and so on\n"
        "@CO\tcitation needed\n"};

    const std::string text = header.ToSam();
    std::shared_ptr<bam_hdr_t> rawData(sam_hdr_parse(text.size(), text.c_str()),
                                       BamHeaderTests::BamHdrDeleter());
    rawData->ignore_sam_err = 0;
    rawData->cigar_tab = nullptr;
    rawData->l_text = text.size();
    rawData->text = static_cast<char*>(calloc(rawData->l_text + 1, 1));
    memcpy(rawData->text, text.c_str(), rawData->l_text);

    const std::string rawText(rawData->text, rawData->l_text);
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
        .PacBioBamVersion("3.0.1")
        .AddReadGroup(rg1)
        .AddReadGroup(rg2)
        .AddReadGroup(rg3)
        .AddSequence(seq1)
        .AddSequence(seq2)
        .AddProgram(prog1)
        .AddComment("ipsum and so on")
        .AddComment("citation needed");

    const std::string expectedText{
        "@HD\tVN:1.1\tSO:queryname\tpb:3.0.1\n"
        "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
        "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
        "@RG\tID:rg1\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:control\tPM:SEQUEL\n"
        "@RG\tID:rg2\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\tPM:SEQUEL\n"
        "@RG\tID:rg3\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\tPM:SEQUEL\n"
        "@PG\tID:_foo_\tPN:ide\n"
        "@CO\tipsum and so on\n"
        "@CO\tcitation needed\n"};

    std::string text = header.ToSam();
    std::shared_ptr<bam_hdr_t> rawData(sam_hdr_parse(text.size(), text.c_str()),
                                       BamHeaderTests::BamHdrDeleter());
    rawData->ignore_sam_err = 0;
    rawData->cigar_tab = nullptr;
    rawData->l_text = text.size();
    rawData->text = static_cast<char*>(calloc(rawData->l_text + 1, 1));
    memcpy(rawData->text, text.c_str(), rawData->l_text);

    const BamHeader newHeader = BamHeader(std::string(rawData->text, rawData->l_text));

    EXPECT_EQ(header.Version(), newHeader.Version());
    EXPECT_EQ(header.SortOrder(), newHeader.SortOrder());
    EXPECT_EQ(header.PacBioBamVersion(), newHeader.PacBioBamVersion());

    text = newHeader.ToSam();
    EXPECT_EQ(expectedText, text);
}

TEST(BamHeaderTest, MergeOk)
{
    const std::string hdrText1{
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:a955def6\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
        "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
        "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.3.0.0.140018;FRAMERATEHZ=75.000000\t"
        "PU:m140918_150013_42139_c100697631700000001823144703261565_s1_p0\t"
        "PM:SEQUEL\n"
        "@PG\tID:bam2bam-0.20.0\tPN:bam2bam\tVN:0.20.0\n"
        "@PG\tID:bax2bam-0.0.2\tPN:bax2bam\tVN:0.0.2\n"
        "@CO\tcomment1\n"};

    const std::string hdrText2{
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;"
        "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;"
        "PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;"
        "PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;"
        "BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;"
        "FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\t"
        "PM:SEQUEL\n"
        "@PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0\n"
        "@PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0\n"
        "@PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0\n"
        "@CO\tcomment2\n"};

    const std::string mergedText{
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:a955def6\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
        "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
        "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.3.0.0.140018;FRAMERATEHZ=75.000000\t"
        "PU:m140918_150013_42139_c100697631700000001823144703261565_s1_p0\t"
        "PM:SEQUEL\n"
        "@RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;"
        "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;"
        "PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;"
        "PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;"
        "BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;"
        "FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\t"
        "PM:SEQUEL\n"
        "@PG\tID:bam2bam-0.20.0\tPN:bam2bam\tVN:0.20.0\n"
        "@PG\tID:bax2bam-0.0.2\tPN:bax2bam\tVN:0.0.2\n"
        "@PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0\n"
        "@PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0\n"
        "@PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0\n"
        "@CO\tcomment1\n"
        "@CO\tcomment2\n"};

    {  // operator+

        const BamHeader header1(hdrText1);
        const BamHeader header2(hdrText2);
        const BamHeader merged = header1 + header2;
        EXPECT_EQ(mergedText, merged.ToSam());

        // also make sure inputs not changed
        EXPECT_EQ(hdrText1, header1.ToSam());
        EXPECT_EQ(hdrText2, header2.ToSam());
    }

    {  // operator+=

        BamHeader header1(hdrText1);
        header1 += BamHeader(hdrText2);
        EXPECT_EQ(mergedText, header1.ToSam());
    }
}

TEST(BamHeaderTest, MergeHandlesDuplicateReadGroups)
{
    const std::string hdrText{
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:a955def6\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
        "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
        "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.3.0.0.140018;FRAMERATEHZ=75.000000\t"
        "PU:m140918_150013_42139_c100697631700000001823144703261565_s1_p0\tPM:SEQUEL\n"
        "@PG\tID:bam2bam-0.20.0\tPN:bam2bam\tVN:0.20.0\n"
        "@PG\tID:bax2bam-0.0.2\tPN:bax2bam\tVN:0.0.2\n"};

    // duplicate @RG:IDs handled ok (i.e. not duplicated in output)
    const BamHeader header1(hdrText);
    const BamHeader header2(hdrText);
    const BamHeader merged = header1 + header2;
    EXPECT_EQ(hdrText, merged.ToSam());
}

TEST(BamHeaderTest, MergeCompatibilityOk)
{
    {  // different @HD:VN - this IS allowed (as of SAT-465, pbbam v0.7.2)
        const std::string hdrText1 = {"@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"};
        const std::string hdrText2 = {"@HD\tVN:1.0\tSO:unknown\tpb:3.0.1\n"};
        const BamHeader header1(hdrText1);
        const BamHeader header2(hdrText2);
        EXPECT_NO_THROW(header1 + header2);
    }

    {  // different @HD:SO
        const std::string hdrText1 = {"@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"};
        const std::string hdrText2 = {"@HD\tVN:1.1\tSO:coordinate\tpb:3.0.1\n"};
        const BamHeader header1(hdrText1);
        const BamHeader header2(hdrText2);
        EXPECT_THROW(header1 + header2, std::runtime_error);
    }

    {  // different @HD:pb - this IS allowed (as of SAT-529, pbbam 0.7.4)
        const std::string hdrText1 = {"@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"};
        const std::string hdrText2 = {"@HD\tVN:1.1\tSO:unknown\tpb:3.0.3\n"};
        const BamHeader header1(hdrText1);
        const BamHeader header2(hdrText2);
        EXPECT_NO_THROW(header1 + header2);
    }

    {  // @SQ list clash
        const std::string hdrText1{
            "@HD\tVN:1.1\tSO:coordinate\tpb:3.0.1\n"
            "@SQ\tSN:foo\tLN:42\n"
            "@SQ\tSN:bar\tLN:24\n"};
        const std::string hdrText2{
            "@HD\tVN:1.1\tSO:coordinate\tpb:3.0.1\n"
            "@SQ\tSN:foo\tLN:42\n"
            "@SQ\tSN:baz\tLN:99\n"};
        const BamHeader header1(hdrText1);
        const BamHeader header2(hdrText2);
        EXPECT_THROW(header1 + header2, std::runtime_error);
    }
}
