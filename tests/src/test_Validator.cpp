#include <pbbam/Validator.h>

#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/ReadGroupInfo.h>

#include <pbcopper/data/Cigar.h>

#include "../src/ValidationErrors.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace ValidatorTests {

static BamRecord makeValidMappedRecord()
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
    impl.CigarData(Data::Cigar{"10="});

    TagCollection tags;
    tags["RG"] = std::string{"db972a04"};
    tags["dq"] = std::string{"2222'$22'2"};
    tags["dt"] = std::string{"NNNNAGNNGN"};
    tags["iq"] = std::string{"(+#1'$#*1&"};
    tags["mq"] = std::string{"&1~51*5&~2"};
    tags["sq"] = std::string{"<32<4<<<<3"};
    tags["ip"] = std::vector<uint8_t>{2, 0, 10, 22, 34, 0, 2, 3, 0, 16};
    tags["np"] = int32_t{1};
    tags["qe"] = int32_t{10};
    tags["qs"] = int32_t{0};
    tags["zm"] = int32_t{54130};
    tags["cx"] = int32_t{2};
    tags["AS"] = int32_t{-3020};
    tags["NM"] = int32_t{134};
    tags["rq"] = static_cast<float>(0.854);
    tags["sn"] = std::vector<float>{2.0, 2.0, 2.0, 2.0};
    impl.Tags(tags);

    return impl;
}

static BamRecord makeValidUnmappedRecord()
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
    tags["RG"] = std::string{"b5482b33"};
    tags["dq"] = std::string{"2222222222"};
    tags["dt"] = std::string{"NNNNNNNNNN"};
    tags["iq"] = std::string{",*11111001"};
    tags["mq"] = std::string{"&47088')34"};
    tags["sq"] = std::string{"8<4<:<6<0<"};
    tags["ip"] = std::vector<uint8_t>{255, 9, 20, 43, 38, 12, 9, 30, 39, 22};
    tags["np"] = int32_t{1};
    tags["qe"] = int32_t{10};
    tags["qs"] = int32_t{0};
    tags["zm"] = int32_t{8};
    tags["cx"] = int32_t{2};
    tags["AS"] = int32_t{-3020};
    tags["NM"] = int32_t{134};
    tags["rq"] = static_cast<float>(0.811);
    tags["sn"] = std::vector<float>{2.0, 2.0, 2.0, 2.0};
    impl.Tags(tags);

    return impl;
}

static ReadGroupInfo makeValidReadGroup()
{
    ReadGroupInfo rg("f5b4ffb6");
    rg.MovieName("movie32");
    rg.ReadType("CCS");
    rg.BindingKit("101-789-500");
    rg.SequencingKit("101-789-300");
    rg.BasecallerVersion("5.0");
    rg.FrameRateHz("100");
    rg.Control("TRUE");
    return rg;
}

// valid, 'starter' objects
static const ReadGroupInfo validReadGroup = makeValidReadGroup();
static const BamRecord validMappedRecord = makeValidMappedRecord();
static const BamRecord validUnmappedRecord = makeValidUnmappedRecord();

}  // namespace ValidatorTests

TEST(ValidatorErrorsTest, default_uses_max)
{
    ValidationErrors errors;
    EXPECT_EQ(ValidationErrors::MAX, errors.MaxNumErrors());
}

TEST(ValidatorErrorsTest, zero_max_num_errors_fallsback_to_max)
{
    ValidationErrors errors(0);
    EXPECT_EQ(ValidationErrors::MAX, errors.MaxNumErrors());
}

TEST(ValidatorErrorsTest, can_specify_max_num_errors)
{
    {  // default - use "no max"
        ValidationErrors errors;
        EXPECT_EQ(ValidationErrors::MAX, errors.MaxNumErrors());
    }
    {  // max of zero doesn't make sense... make equivalent to "no max"
        ValidationErrors errors(0);
        EXPECT_EQ(ValidationErrors::MAX, errors.MaxNumErrors());
    }
    {  // max = 1
        ValidationErrors errors(1);
        EXPECT_EQ(1, errors.MaxNumErrors());
    }
    {  // max = 10
        ValidationErrors errors(10);
        EXPECT_EQ(10, errors.MaxNumErrors());
    }
}

TEST(ValidatorErrorsTest, throws_when_max_errors_reached)
{
    {
        ValidationErrors errors(1);
        EXPECT_THROW(errors.AddFileError("foo", "you"), ValidationException);
    }
    {
        ValidationErrors errors(2);
        errors.AddFileError("foo", "you");
        EXPECT_THROW(errors.AddFileError("foo", "me"), ValidationException);
    }
}

TEST(ValidatorErrorsTest, constructs_expected_validation_exception)
{
    const std::string error1 = "error1";
    const std::string error2 = "error2";

    try {

        ValidationErrors errors(4);
        errors.AddFileError("path/to/foo.bam", error1);
        errors.AddFileError("path/to/foo.bam", error2);
        errors.AddReadGroupError("deadbeef", "invalid sequencing chemistry combination detected");
        errors.AddRecordError(
            "m140906_231018_42161_c100676332550000001823129611271486_s1_p0/8/0_10",
            "MergeQV does not match expected length");

    } catch (ValidationException& e) {

        EXPECT_EQ(1, e.FileErrors().size());                        // only 1 file
        EXPECT_EQ(2, e.FileErrors().at("path/to/foo.bam").size());  // 2 errors for this file
        EXPECT_EQ(1, e.ReadGroupErrors().size());
        EXPECT_EQ(1, e.RecordErrors().size());
    }
}

TEST(BAM_Validator, success_on_valid_read_group)
{
    ASSERT_NO_THROW(Validator::Validate(ValidatorTests::validReadGroup));
}

TEST(BAM_Validator, reports_missing_read_group_components)
{
    {  // missing ID
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.Id("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // missing movie name
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.MovieName("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // missing read type
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.ReadType("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // missing binding kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BindingKit("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // missing sequencing kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.SequencingKit("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // missing basecaller version
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BasecallerVersion("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // missing frame rate
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.FrameRateHz("");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
}

TEST(BAM_Validator, reports_invalid_read_group_data)
{
    {  // mismatch expected ID vs stored ID - change ID
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.Id("deadbeef");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // mismatch expected ID vs stored ID - change read type
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.ReadType("SUBREAD");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // mismatch expected ID vs stored ID - change movie name
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.MovieName("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // unknown read type
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.ReadType("FOO");

        // recompute ID so we're only checking the new read type, not read ID
        rg.Id(MakeReadGroupId(rg.MovieName(), rg.ReadType()));

        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // invalid chemistry triple - change binding kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BindingKit("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // invalid chemistry triple - change sequencing kit
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.SequencingKit("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // invalid chemistry triple - change basecaller version
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.BasecallerVersion("0.42");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
    {  // non-numeric frame rate
        ReadGroupInfo rg = ValidatorTests::validReadGroup;
        rg.FrameRateHz("foo");
        EXPECT_THROW(Validator::Validate(rg), ValidationException);
        EXPECT_FALSE(Validator::IsValid(rg));
    }
}

TEST(BAM_Validator, success_on_valid_header)
{
    static const BamHeader validMappedHeader{
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.7\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:"
        "734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000"
        "\tPU:m64004_190414_193017\tPM:SEQUELII\n"};

    static const BamHeader validUnmappedHeader{
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000\t"
        "PU:m64004_190414_193017\tPM:SEQUELII\n"};

    ASSERT_NO_THROW(Validator::Validate(validMappedHeader));
    ASSERT_NO_THROW(Validator::Validate(validUnmappedHeader));
}

TEST(BAM_Validator, success_on_valid_revio_header)
{
    static const BamHeader validMappedHeader{
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.7\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:"
        "734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000"
        "\tPU:m64004_190414_193017\tPM:REVIO\n"};

    static const BamHeader validUnmappedHeader{
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000\t"
        "PU:m64004_190414_193017\tPM:REVIO\n"};

    ASSERT_NO_THROW(Validator::Validate(validMappedHeader));
    ASSERT_NO_THROW(Validator::Validate(validUnmappedHeader));
}

TEST(BAM_Validator, reports_invalid_header_data)
{
    static const BamHeader validMappedHeader{
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.7\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:"
        "734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000"
        "\tPU:m64004_190414_193017\tPM:SEQUELII\n"};

    {  // invalid SAM version - non-numeric
        BamHeader header = validMappedHeader.DeepCopy();
        header.Version("foo");
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {  // invalid SAM version - negative version numbers
        BamHeader header = validMappedHeader.DeepCopy();
        header.Version("-1.4.0");
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {  // invalid sort order
        BamHeader header = validMappedHeader.DeepCopy();
        header.SortOrder("not_a_valid_sort_order");
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }

    // invalid PacBioBamVersion numbers (non-numeric, negative, earlier than min)
    // already throw when you try to set them... so we have to catch & ignore
    // initial exception to get to validator

    {  // invalid PacBioBAM version - non-numeric
        BamHeader header = validMappedHeader.DeepCopy();
        try {
            header.PacBioBamVersion("foo");
        } catch (...) {
        }
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {  // invalid PacBioBAM version - negative version numbers
        BamHeader header = validMappedHeader.DeepCopy();
        try {
            header.PacBioBamVersion("-1.4.0");
        } catch (...) {
        }
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
    {  // invalid PacBioBAM version - earlier than minimum allowed
        BamHeader header = validMappedHeader.DeepCopy();
        try {
            header.PacBioBamVersion("3.0.0");
        } catch (...) {
        }
        EXPECT_THROW(Validator::Validate(header), ValidationException);
        EXPECT_FALSE(Validator::IsValid(header));
    }
}

TEST(BAM_Validator, success_on_valid_bam_record)
{
    static const BamHeader validMappedHeader{
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.7\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:"
        "734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000"
        "\tPU:m64004_190414_193017\tPM:SEQUELII\n"};
    BamRecord record(ValidatorTests::validMappedRecord);
    record.header_ = validMappedHeader;
    ASSERT_NO_THROW(Validator::Validate(record));
}

static void ModifyTag(BamRecord* record, const std::string& tagName, const Tag& tag)
{
    if (record->Impl().HasTag(tagName)) {
        record->Impl().EditTag(tagName, tag);
    } else {
        record->Impl().AddTag(tagName, tag);
    }
}

static void CheckInvalidTagLength(const std::string& tagName, const Tag& tag)
{
    static const BamHeader validUnmappedHeader{
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000\t"
        "PU:m64004_190414_193017\tPM:SEQUELII\n"};
    BamRecord record(ValidatorTests::validUnmappedRecord);
    record.header_ = validUnmappedHeader;

    ModifyTag(&record, tagName, tag);

    EXPECT_THROW(Validator::Validate(record), ValidationException);
    EXPECT_FALSE(Validator::IsValid(record));
}

TEST(BAM_Validator, reports_invalid_tag_lengths)
{
    static const BamHeader validUnmappedHeader{
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000\t"
        "PU:m64004_190414_193017\tPM:SEQUELII\n"};

    // make these "variable-length" SEQ/tags too short for the read's stated
    // queryStart/queryEnd

    {  // SEQ
        BamRecord record(ValidatorTests::validUnmappedRecord);
        record.header_ = validUnmappedHeader;
        record.Impl().SetSequenceAndQualities("AA");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }

    CheckInvalidTagLength("dq", Data::QualityValues("@@").Fastq());  // DeletionQV
    CheckInvalidTagLength("iq", Data::QualityValues("@@").Fastq());  // InsertionQV
    CheckInvalidTagLength("mq", Data::QualityValues("@@").Fastq());  // MergeQV
    CheckInvalidTagLength("sq", Data::QualityValues("@@").Fastq());  // SubstitutionQV
    CheckInvalidTagLength("dt", std::string("AA"));                  // DeletionTag
    CheckInvalidTagLength("st", std::string("AA"));                  // SubstitutionTag

    const Data::Frames f{{42, 42, 42}};
    const auto& frames = f.Data();
    CheckInvalidTagLength("ip", frames);  // IPD

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

TEST(BAM_Validator, reports_invalid_tag_data)
{
    static const BamHeader validMappedHeader{
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.7\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:"
        "734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000"
        "\tPU:m64004_190414_193017\tPM:SEQUELII\n"};

    {  // missing qe
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("qe");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {  // missing qs
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("qs");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {  // queryStart should be < queryEnd
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
    {  // missing np
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("np");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {  // numPasses for SUBREAD type records should be 1
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.NumPasses(42);
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {  // missing sn
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().RemoveTag("sn");
        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
}

TEST(BAM_Validator, reports_invalid_mapped_record_data)
{
    static const BamHeader validMappedHeader{
        "@HD\tVN:1.5\tSO:coordinate\tpb:3.0.7\n"
        "@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:"
        "734d5f3b2859595f4bd87a2fe6b7389b\n"
        "@RG\tID:db972a04\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000"
        "\tPU:m64004_190414_193017\tPM:SEQUELII\n"};

    {  // mapped record should have valid refID
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().ReferenceId(-1);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {  // mapped record should have valid position
        BamRecord record(ValidatorTests::validMappedRecord);
        record.header_ = validMappedHeader;
        record.Impl().Position(-1);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
}

TEST(BAM_Validator, reports_invalid_unmapped_record_data)
{
    static const BamHeader validUnmappedHeader{
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n"
        "@RG\tID:b5482b33\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
        "BINDINGKIT=101-717-300;SEQUENCINGKIT=101-644-500;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100."
        "000000\t"
        "PU:m64004_190414_193017\tPM:SEQUELII\n"};

    {  // unmapped should have no refID
        BamRecord record(ValidatorTests::validUnmappedRecord);
        record.header_ = validUnmappedHeader;
        record.Impl().ReferenceId(0);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
    {  // unmapped should have no position
        BamRecord record(ValidatorTests::validUnmappedRecord);
        record.header_ = validUnmappedHeader;
        record.Impl().Position(42);

        EXPECT_THROW(Validator::Validate(record), ValidationException);
        EXPECT_FALSE(Validator::IsValid(record));
    }
}
