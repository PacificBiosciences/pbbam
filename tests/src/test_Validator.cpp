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

#include <gtest/gtest.h>

#ifdef PBBAM_TESTING
#define private public
#endif

#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/Cigar.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/Validator.h>

#include "../src/StringUtils.h"
#include "../src/ValidationErrors.h"

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace ValidatorTests {

static BamRecord makeValidMappedRecord(void)
{
    BamRecordImpl impl;
    impl.Bin(4680);
    impl.Flag(2);
    impl.InsertSize(0);
    impl.MapQuality(10);
    impl.MatePosition(-1);
    impl.MateReferenceId(-1);
    impl.Name("movie1/54130/0_10");
    impl.Position(1);
    impl.ReferenceId(0);
    impl.SetMapped(true);
    impl.SetSequenceAndQualities("AATGAGGAGA");
    impl.CigarData(Cigar{ "10=" });

    TagCollection tags;
    tags["RG"] = string{ "3f58e5b8" };
    tags["dq"] = string{ "2222'$22'2" };
    tags["dt"] = string{ "NNNNAGNNGN" };
    tags["iq"] = string{ "(+#1'$#*1&" };
    tags["mq"] = string{ "&1~51*5&~2" };
    tags["sq"] = string{ "<32<4<<<<3" };
    tags["ip"] = vector<uint8_t>{ 2,0,10,22,34,0,2,3,0,16 };
    tags["np"] = static_cast<int32_t>(1);
    tags["qe"] = static_cast<int32_t>(10);
    tags["qs"] = static_cast<int32_t>(0);
    tags["zm"] = static_cast<int32_t>(54130);
    tags["cx"] = static_cast<int32_t>(2);
    tags["AS"] = static_cast<int32_t>(-3020);
    tags["NM"] = static_cast<int32_t>(134);
    tags["rq"] = static_cast<float>(0.854);
    tags["sn"] = vector<float>{ 2.0,2.0,2.0,2.0 };
    impl.Tags(tags);

    return BamRecord(impl);
}

static BamRecord makeValidUnmappedRecord(void)
{
    BamRecordImpl impl;
    impl.Bin(4680);
    impl.Flag(4);
    impl.InsertSize(0);
    impl.MapQuality(10);
    impl.MatePosition(-1);
    impl.MateReferenceId(-1);
    impl.Name("m140906_231018_42161_c100676332550000001823129611271486_s1_p0/8/0_10");
    impl.Position(-1);
    impl.ReferenceId(-1);
    impl.SetSequenceAndQualities("AATGAGGAGA");

    TagCollection tags;
    tags["RG"] = string{ "b5482b33" };
    tags["dq"] = string{ "2222222222" };
    tags["dt"] = string{ "NNNNNNNNNN" };
    tags["iq"] = string{ ",*11111001" };
    tags["mq"] = string{ "&47088')34" };
    tags["sq"] = string{ "8<4<:<6<0<" };
    tags["ip"] = vector<uint8_t>{ 255,9,20,43,38,12,9,30,39,22 };
    tags["np"] = static_cast<int32_t>(1);
    tags["qe"] = static_cast<int32_t>(10);
    tags["qs"] = static_cast<int32_t>(0);
    tags["zm"] = static_cast<int32_t>(8);
    tags["cx"] = static_cast<int32_t>(2);
    tags["AS"] = static_cast<int32_t>(-3020);
    tags["NM"] = static_cast<int32_t>(134);
    tags["rq"] = static_cast<float>(0.811);
    tags["sn"] = vector<float>{ 2.0,2.0,2.0,2.0 };
    impl.Tags(tags);

    return BamRecord(impl);
}

static ReadGroupInfo makeValidReadGroup(void)
{
    ReadGroupInfo rg("f5b4ffb6");
    rg.MovieName("movie32");
    rg.ReadType("CCS");
    rg.BindingKit("100372700");
    rg.SequencingKit("100612400");
    rg.BasecallerVersion("2.3");
    rg.FrameRateHz("100");
    rg.Control("TRUE");
    return rg;
}

// valid, 'starter' objects
static const ReadGroupInfo validReadGroup       = makeValidReadGroup();
static const BamRecord     validMappedRecord    = makeValidMappedRecord();
static const BamRecord     validUnmappedRecord  = makeValidUnmappedRecord();

} // namespace ValidatorTests

TEST(ValidatorErrorsTest, SetMaxNumErrors)
{
    {   // default - use "no max"
        internal::ValidationErrors errors;
        EXPECT_EQ(internal::ValidationErrors::MAX, errors.maxNumErrors_);
    }
    {   // max of zero doesn't make sense... make equivalent to "no max"
        internal::ValidationErrors errors(0);
        EXPECT_EQ(internal::ValidationErrors::MAX, errors.maxNumErrors_);
    }
    {   // max = 1
        internal::ValidationErrors errors(1);
        EXPECT_EQ(1, errors.maxNumErrors_);
    }
    {   // max = 10
        internal::ValidationErrors errors(10);
        EXPECT_EQ(10, errors.maxNumErrors_);
    }
}

TEST(ValidatorErrorsTest, ThrowOnMaxReached)
{
    {
        internal::ValidationErrors errors(1);
        EXPECT_THROW(errors.AddFileError("foo", "you"), ValidationException);
    }
    {
        internal::ValidationErrors errors(2);
        errors.AddFileError("foo", "you");
        EXPECT_THROW(errors.AddFileError("foo", "me"), ValidationException);
    }
}

TEST(ValidatorErrorsTest, ExceptionFromResults)
{
    const string error1 = "error1";
    const string error2 = "error2";

    try {

        internal::ValidationErrors errors(4);
        errors.AddFileError("path/to/foo.bam", error1);
        errors.AddFileError("path/to/foo.bam", error2);
        errors.AddReadGroupError("deadbeef", "invalid sequencing chemistry combination detected");
        errors.AddRecordError("m140906_231018_42161_c100676332550000001823129611271486_s1_p0/8/0_10",
                              "MergeQV does not match expected length");

    } catch (ValidationException& e) {

        EXPECT_EQ(1, e.FileErrors().size());                       // only 1 file
        EXPECT_EQ(2, e.FileErrors().at("path/to/foo.bam").size()); // 2 errors for this file
        EXPECT_EQ(1, e.ReadGroupErrors().size());
        EXPECT_EQ(1, e.RecordErrors().size());
    }
}

TEST(ValidatorTest, ValidReadGroup)
{
    ASSERT_NO_THROW(Validator::Validate(ValidatorTests::validReadGroup));
}

TEST(ValidatorTest, ReadGroupRequiredComponents)
{
    {   // missing ID
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.Id("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // missing movie name
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.MovieName("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // missing read type
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.ReadType("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // missing binding kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BindingKit("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // missing sequencing kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.SequencingKit("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // missing basecaller version
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BasecallerVersion("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // missing frame rate
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.FrameRateHz("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
}

TEST(ValidatorTest, ReadGroupValues)
{
    {   // mismatch expected ID vs stored ID - change ID
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.Id("deadbeef");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // mismatch expected ID vs stored ID - change read type
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.ReadType("SUBREAD");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // mismatch expected ID vs stored ID - change movie name
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.MovieName("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // unknown read type
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.ReadType("FOO");

        // recompute ID so we're only checking the new read type, not read ID
        rg.Id( MakeReadGroupId(rg.MovieName(), rg.ReadType()) );

        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // invalid chemistry triple - change binding kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BindingKit("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // invalid chemistry triple - change sequencing kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.SequencingKit("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {   // invalid chemistry triple - change basecaller version
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BasecallerVersion("0.42");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    { // non-numeric frame rate
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.FrameRateHz("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
}

TEST(ValidatorTest, ValidHeader)
{
    const BamHeader validMappedHeader {
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.1\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:3f58e5b8\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;"
              "MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BASECALLERVERSION=2.1;"
              "FRAMERATEHZ=75.000000;BINDINGKIT=100356300;SEQUENCINGKIT=100356200"
              "\tPU:movie1\n"
    };

    const BamHeader validUnmappedHeader {
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
            "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
            "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.1;FRAMERATEHZ=75.000000\t"
            "PU:m140906_231018_42161_c100676332550000001823129611271486_s1_p0\n"
    };

    ASSERT_NO_THROW(Validator::Validate(validMappedHeader));
    ASSERT_NO_THROW(Validator::Validate(validUnmappedHeader));
}

TEST(ValidatorTest, ValidateHeader)
{
    const BamHeader validMappedHeader {
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.1\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:3f58e5b8\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;"
              "MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BASECALLERVERSION=2.1;"
              "FRAMERATEHZ=75.000000;BINDINGKIT=100356300;SEQUENCINGKIT=100356200"
              "\tPU:movie1\n"
    };

    {   // invalid SAM version - non-numeric
        BamHeader header = validMappedHeader.DeepCopy();
        header.Version("foo");
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {   // invalid SAM version - negative version numbers
        BamHeader header = validMappedHeader.DeepCopy();
        header.Version("-1.4.0");
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {   // invalid sort order
        BamHeader header = validMappedHeader.DeepCopy();
        header.SortOrder("not_a_valid_sort_order");
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }

    // invalid PacBioBamVersion numbers (non-numeric, negative, earlier than min)
    // already throw when you try to set them... so we have to catch & ignore
    // initial exception to get to validator

    {   // invalid PacBioBAM version - non-numeric
        BamHeader header = validMappedHeader.DeepCopy();
        try {
            header.PacBioBamVersion("foo");
        } catch (...) { }
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {   // invalid PacBioBAM version - negative version numbers
        BamHeader header = validMappedHeader.DeepCopy();
        try {
            header.PacBioBamVersion("-1.4.0");
        } catch (...) { }
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {   // invalid PacBioBAM version - earlier than minimum allowed
        BamHeader header = validMappedHeader.DeepCopy();
        try {
            header.PacBioBamVersion("3.0.0");
        } catch (...) { }
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
}

TEST(ValidatorTest, ValidRecord)
{
    const BamHeader validMappedHeader {
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.1\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:3f58e5b8\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;"
              "MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BASECALLERVERSION=2.1;"
              "FRAMERATEHZ=75.000000;BINDINGKIT=100356300;SEQUENCINGKIT=100356200"
              "\tPU:movie1\n"
    };
    BamRecord record(ValidatorTests::validMappedRecord);
    record.header_ = validMappedHeader;
    ASSERT_NO_THROW(Validator::Validate(record));
}

static inline
void ModifyTag(BamRecord* record,
               const std::string& tagName,
               const Tag& tag)
{
    if (record->Impl().HasTag(tagName))
        record->Impl().EditTag(tagName, tag);
    else
        record->Impl().AddTag(tagName, tag);
}

static inline
void CheckInvalidTagLength(const std::string& tagName, const Tag& tag)
{
    static const BamHeader validUnmappedHeader {
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
            "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
            "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.1;FRAMERATEHZ=75.000000\t"
            "PU:m140906_231018_42161_c100676332550000001823129611271486_s1_p0\n"
    };
    BamRecord record(ValidatorTests::validUnmappedRecord);
    record.header_ = validUnmappedHeader;

    ModifyTag(&record, tagName, tag);

    EXPECT_THROW(Validator::Validate(record), ValidationException);
    EXPECT_FALSE(Validator::IsValid(record));
}

TEST(ValidatorTest, TagDataLengths)
{
    const BamHeader validUnmappedHeader {
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
            "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
            "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.1;FRAMERATEHZ=75.000000\t"
            "PU:m140906_231018_42161_c100676332550000001823129611271486_s1_p0\n"
    };

    // make these "variable-length" SEQ/tags too short for the read's stated
    // queryStart/queryEnd

    {   // SEQ
        BamRecord record(ValidatorTests::validUnmappedRecord);
        record.header_ = validUnmappedHeader;
        record.Impl().SetSequenceAndQualities("AA");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }

    CheckInvalidTagLength("dq", QualityValues("@@").Fastq());  // DeletionQV
    CheckInvalidTagLength("iq", QualityValues("@@").Fastq());  // InsertionQV
    CheckInvalidTagLength("mq", QualityValues("@@").Fastq());  // MergeQV
    CheckInvalidTagLength("sq", QualityValues("@@").Fastq());  // SubstitutionQV
    CheckInvalidTagLength("dt", string("AA")); // DeletionTag
    CheckInvalidTagLength("st", string("AA")); // SubstitutionTag

    const auto& f = Frames{ vector<uint16_t>{42, 42, 42} };
    const auto& frames = f.Data();
    CheckInvalidTagLength("ip", frames); // IPD

    // NOTE: disabling "internal" tag checks for now, only checking "standard"
    //       PacBioBAM tags

//    const auto& pulses = vector<uint16_t>{42, 42, 42};
//    CheckInvalidTagLength("pv", QualityValues("@@").Fastq());  // AltLabelQV
//    CheckInvalidTagLength("pq", QualityValues("@@").Fastq());  // LabelQV
//    CheckInvalidTagLength("pg", QualityValues("@@").Fastq());  // PulseMergeQv
//    CheckInvalidTagLength("pt", string("AA")); // AltLabelTag
//    CheckInvalidTagLength("pc", string("AA")); // PulseCall
//    CheckInvalidTagLength("pd", frames); // PrePulseFrames
//    CheckInvalidTagLength("px", frames); // PulseCallWidth
//    CheckInvalidTagLength("pw", frames); // PulseWidth
//    CheckInvalidTagLength("pa", pulses); // Pkmean
//    CheckInvalidTagLength("ps", pulses); // Pkmean2
//    CheckInvalidTagLength("pm", pulses); // Pkmid
//    CheckInvalidTagLength("pi", pulses); // Pkmid2
}

TEST(ValidatorTest, TagDataValues)
{
    const BamHeader validMappedHeader {
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.1\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:3f58e5b8\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;"
              "MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BASECALLERVERSION=2.1;"
              "FRAMERATEHZ=75.000000;BINDINGKIT=100356300;SEQUENCINGKIT=100356200"
              "\tPU:movie1\n"
    };

    {   // missing qe
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("qe");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {   // missing qs
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("qs");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {   // queryStart should be < queryEnd
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.QueryStart(10);
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {  // missing zm
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("zm");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {   // missing np
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("np");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {   // numPasses for SUBREAD type records should be 1
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.NumPasses(42);
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {   // missing sn
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("sn");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
}

TEST(ValidatorTest, MappedRecords)
{
    const BamHeader validMappedHeader {
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.1\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
            "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
            "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.1;FRAMERATEHZ=75.000000\t"
            "PU:m140906_231018_42161_c100676332550000001823129611271486_s1_p0\n"
    };

    {   // mapped record should have valid refID
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().ReferenceId(-1);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {   // mapped record should have valid position
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().Position(-1);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }

}

TEST(ValidatorTest, UnmappedRecords)
{
    const BamHeader validUnmappedHeader {
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;"
            "InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;BINDINGKIT=100356300;"
            "SEQUENCINGKIT=100356200;BASECALLERVERSION=2.1;FRAMERATEHZ=75.000000\t"
            "PU:m140906_231018_42161_c100676332550000001823129611271486_s1_p0\n"
    };

    {   // unmapped should have no refID
        BamRecord record(ValidatorTests::validUnmappedRecord);
        record.header_ = validUnmappedHeader;
        record.Impl().ReferenceId(0);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {   // unmapped should have no position
        BamRecord record(ValidatorTests::validUnmappedRecord);
        record.header_ = validUnmappedHeader;
        record.Impl().Position(42);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
}
