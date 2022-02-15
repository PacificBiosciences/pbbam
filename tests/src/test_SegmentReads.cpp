#include "PbbamTestData.h"

#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/ZmwGroupQuery.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include <cstdio>

using namespace PacBio;

namespace SegmentReadTests {

const std::string BasicSegmentBamFn{BAM::PbbamTestsConfig::Data_Dir +
                                    "/segment/basic.ccs.segments.bam"};

}  // namespace SegmentReadTests

TEST(BAM_SegmentReads, can_query_header_read_groups)
{
    const BAM::BamFile file{SegmentReadTests::BasicSegmentBamFn};

    const BAM::BamHeader& header = file.Header();
    EXPECT_TRUE(header.HasReadGroup("e51ee4ef"));
    const BAM::ReadGroupInfo& rg = header.ReadGroup("e51ee4ef");
    EXPECT_TRUE(rg.IsSegment());
    EXPECT_EQ(rg.SegmentSource().value(), "CCS");

    BAM::BamRecord record{header};
    record.ReadGroup(rg);
    EXPECT_TRUE(record.IsSegment());
    EXPECT_TRUE(record.Type() == BAM::RecordType::SEGMENT);
}

TEST(BAM_SegmentReads, can_query_edit_read_group_objects)
{
    const std::string ccsRgSam{
        "@RG\tID:9e129d4c\tPL:PACBIO\tDS:READTYPE=CCS;"
        "BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;"
        "FRAMERATEHZ=100.000000\tLB:AML_Mas-seq_SD\tPU:m64013e_211031_055434\t"
        "SM:AML_Mas-seq_SD\tPM:SEQUELII\tCM:S/P5-C2/5.0-8M"};
    const std::string segmentRgSam{
        "@RG\tID:e51ee4ef\tPL:PACBIO\tDS:READTYPE=SEGMENT;SOURCE=CCS;"
        "BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;"
        "FRAMERATEHZ=100.000000\tLB:AML_Mas-seq_SD\tPU:m64013e_211031_055434\t"
        "SM:AML_Mas-seq_SD\tPM:SEQUELII\tCM:S/P5-C2/5.0-8M"};

    BAM::ReadGroupInfo rg = BAM::ReadGroupInfo::FromSam(ccsRgSam);
    EXPECT_EQ(rg.ReadType(), "CCS");

    rg.MakeSegment();
    EXPECT_EQ(rg.ReadType(), "SEGMENT");
    EXPECT_TRUE(rg.IsSegment());
    EXPECT_EQ(rg.SegmentSource().value(), "CCS");
    EXPECT_EQ(rg.ToSam(), segmentRgSam);

    rg.RevertSegment();
    EXPECT_EQ(rg.ReadType(), "CCS");
    EXPECT_FALSE(rg.IsSegment());
    EXPECT_FALSE(rg.SegmentSource().has_value());
    EXPECT_EQ(rg.ToSam(), ccsRgSam);
}

TEST(BAM_SegmentReads, can_query_edit_segment_read_tags)
{
    BAM::BamRecord b;

    EXPECT_FALSE(b.HasSegmentIndex());
    EXPECT_FALSE(b.HasSegmentLeftAdapterIndex());
    EXPECT_FALSE(b.HasSegmentRightAdapterIndex());
    EXPECT_FALSE(b.HasSegmentSupplementalData());
    EXPECT_THROW(b.SegmentIndex(), std::runtime_error);
    EXPECT_THROW(b.SegmentLeftAdapterIndex(), std::runtime_error);
    EXPECT_THROW(b.SegmentRightAdapterIndex(), std::runtime_error);
    EXPECT_THROW(b.SegmentSupplementalData(), std::runtime_error);

    const int segmentIndex = 8;
    const int segmentLeftAdapterIndex = 1;
    const int segmentRightAdapterIndex = 3;
    b.SegmentIndex(segmentIndex)
        .SegmentLeftAdapterIndex(segmentLeftAdapterIndex)
        .SegmentRightAdapterIndex(segmentRightAdapterIndex);

    JSON::Json segmentSupplementalJson;
    segmentSupplementalJson["left"] = JSON::Json::object();
    segmentSupplementalJson["left"]["sequence"] = "ACCCGATCA";
    segmentSupplementalJson["left"]["class"] = "RANDOM";
    segmentSupplementalJson["left"]["adapter"] = "";
    segmentSupplementalJson["right"] = JSON::Json::object();
    segmentSupplementalJson["right"]["sequence"] = "GGTTAATTA";
    segmentSupplementalJson["right"]["class"] = "FAILED";
    segmentSupplementalJson["right"]["adapter"] = "ACCCGTAA";
    b.SegmentSupplementalData(segmentSupplementalJson);

    ASSERT_TRUE(b.HasSegmentIndex());
    ASSERT_TRUE(b.HasSegmentLeftAdapterIndex());
    ASSERT_TRUE(b.HasSegmentRightAdapterIndex());
    ASSERT_TRUE(b.HasSegmentSupplementalData());
    EXPECT_EQ(segmentIndex, b.SegmentIndex());
    EXPECT_EQ(segmentLeftAdapterIndex, b.SegmentLeftAdapterIndex());
    EXPECT_EQ(segmentRightAdapterIndex, b.SegmentRightAdapterIndex());
    EXPECT_EQ(segmentSupplementalJson, b.SegmentSupplementalData());

    // raw storage is binary
    const auto supplementalTag = b.Impl().TagValue("ds");
    EXPECT_TRUE(supplementalTag.IsUInt8Array());
}

// m64013e_211031_055434/9830824/ccs/15_1545
// m64013e_211031_055434/9830824/ccs/1561_3070
// m64013e_211031_055434/9830824/ccs/3086_3831
// m64013e_211031_055434/9830824/ccs/3847_5145
// m64013e_211031_055434/9830824/ccs/5161_5980
// m64013e_211031_055434/9830824/ccs/5996_6678
// m64013e_211031_055434/9830824/ccs/6694_7324
// m64013e_211031_055434/9830824/ccs/7340_8385
// m64013e_211031_055434/9830824/ccs/8401_11035
// m64013e_211031_055434/9830824/ccs/11051_11761
// m64013e_211031_055434/9830824/ccs/11777_12282
// m64013e_211031_055434/9830824/ccs/12298_13333
// m64013e_211031_055434/9830824/ccs/13349_14428
// m64013e_211031_055434/9830824/ccs/14444_15112
// m64013e_211031_055434/9830824/ccs/15128_16105
// m64013e_211031_055434/4280389/ccs/15_1545
// m64013e_211031_055434/4280389/ccs/1561_3070
// m64013e_211031_055434/4280389/ccs/3086_3831
// m64013e_211031_055434/4280389/ccs/3847_5145
// m64013e_211031_055434/4280389/ccs/5161_5980
// m64013e_211031_055434/4280389/ccs/5996_6678

TEST(BAM_SegmentReads, pbi_filter_query_can_filter_ccs_segment_records_by_qname)
{
    const std::vector<std::string> qnames{"m64013e_211031_055434/9830824/ccs/7340_8385",
                                          "m64013e_211031_055434/4280389/ccs/3847_5145"};

    std::vector<BAM::BamRecord> records;
    BAM::PbiFilterQuery query{BAM::PbiQueryNameFilter{qnames}, SegmentReadTests::BasicSegmentBamFn};
    for (const auto& record : query) {
        records.push_back(record);
    }
    ASSERT_EQ(records.size(), 2);

    EXPECT_TRUE(std::all_of(records.cbegin(), records.cend(), [](const auto& record) {
        return record.IsSegment() &&
               ((record.FullName() == "m64013e_211031_055434/9830824/ccs/7340_8385") ||
                (record.FullName() == "m64013e_211031_055434/4280389/ccs/3847_5145"));
    }));
}

TEST(BAM_SegmentReads, pbi_filter_query_can_filter_ccs_segment_records_by_zmw)
{
    std::vector<BAM::BamRecord> records;
    BAM::PbiFilterQuery query{BAM::PbiZmwFilter{4280389}, SegmentReadTests::BasicSegmentBamFn};
    for (const auto& record : query) {
        records.push_back(record);
    }
    ASSERT_EQ(records.size(), 6);

    EXPECT_TRUE(std::all_of(records.cbegin(), records.cend(), [](const auto& record) {
        return record.IsSegment() && (record.HoleNumber() == 4280389);
    }));
}

TEST(BAM_SegmentReads, can_get_segments_with_zmw_group_query)
{
    std::vector<std::vector<BAM::BamRecord>> zmws;

    BAM::ZmwGroupQuery query{SegmentReadTests::BasicSegmentBamFn};
    for (const auto& zmw : query) {
        zmws.push_back(zmw);
    }
    ASSERT_EQ(zmws.size(), 2);

    ASSERT_EQ(zmws.at(0).size(), 15);
    EXPECT_TRUE(std::all_of(zmws.at(0).cbegin(), zmws.at(0).cend(), [](const auto& record) {
        return record.IsSegment() && (record.HoleNumber() == 9830824);
    }));

    ASSERT_EQ(zmws.at(1).size(), 6);
    EXPECT_TRUE(std::all_of(zmws.at(1).cbegin(), zmws.at(1).cend(), [](const auto& record) {
        return record.IsSegment() && (record.HoleNumber() == 4280389);
    }));
}

TEST(BAM_SegmentReads, can_get_filtered_segments_with_zmw_group_query)
{
    std::vector<std::vector<BAM::BamRecord>> zmws;

    BAM::ZmwGroupQuery query{{4280389}, SegmentReadTests::BasicSegmentBamFn};
    for (const auto& zmw : query) {
        zmws.push_back(zmw);
    }
    ASSERT_EQ(zmws.size(), 1);

    ASSERT_EQ(zmws.at(0).size(), 6);
    EXPECT_TRUE(std::all_of(zmws.at(0).cbegin(), zmws.at(0).cend(), [](const auto& record) {
        return record.IsSegment() && (record.HoleNumber() == 4280389);
    }));
}

TEST(BAM_SegmentReads, can_make_and_revert_segment_records)
{
    const std::string tempCcsBamFn{BAM::PbbamTestsConfig::GeneratedData_Dir +
                                   "/segment_test.ccs.bam"};
    const std::string tempSegmentBamFn{BAM::PbbamTestsConfig::GeneratedData_Dir +
                                       "/segment_test.segments.bam"};

    // NOTE: This test will not do any ZMW-level work (e.g. stitching), just
    //       making sure BAM record modifications carry across conversions:
    //
    //       SEGMENT -> CCS -> SEGMENT
    //

    // revert to CCS
    {
        // check initial segment header, then revert to CCS read groups
        BAM::BamFile file{SegmentReadTests::BasicSegmentBamFn};
        const BAM::BamHeader& originalHeader = file.Header();
        auto readGroups = originalHeader.ReadGroups();
        ASSERT_EQ(readGroups.size(), 1);
        EXPECT_EQ(readGroups[0].Id(), "e51ee4ef");
        EXPECT_TRUE(readGroups[0].IsSegment());

        BAM::BamHeader ccsHeader = file.Header().DeepCopy();
        readGroups[0].RevertSegment();
        EXPECT_FALSE(readGroups[0].IsSegment());
        EXPECT_EQ(readGroups[0].Id(), "9e129d4c");
        ccsHeader.ReadGroups(readGroups);

        // update records and write to new file
        BAM::BamWriter ccsWriter{tempCcsBamFn, ccsHeader};
        BAM::EntireFileQuery query{SegmentReadTests::BasicSegmentBamFn};
        for (auto& record : query) {
            EXPECT_TRUE(record.IsSegment());
            record.header_ = ccsHeader;
            record.ReadGroup(readGroups[0]);
            EXPECT_FALSE(record.IsSegment());
            ccsWriter.Write(record);
        }
    }
    {
        // check new CCS BAM header, then make segment read groups
        BAM::BamFile file{tempCcsBamFn};
        const BAM::BamHeader& originalHeader = file.Header();
        auto readGroups = originalHeader.ReadGroups();
        ASSERT_EQ(readGroups.size(), 1);
        EXPECT_EQ(readGroups[0].Id(), "9e129d4c");
        EXPECT_FALSE(readGroups[0].IsSegment());

        BAM::BamHeader segmentHeader = file.Header().DeepCopy();
        readGroups[0].MakeSegment();
        EXPECT_TRUE(readGroups[0].IsSegment());
        EXPECT_EQ(readGroups[0].Id(), "e51ee4ef");
        segmentHeader.ReadGroups(readGroups);

        // update records and write to new file
        BAM::BamWriter segmentWriter{tempSegmentBamFn, segmentHeader};
        BAM::EntireFileQuery query{tempCcsBamFn};
        for (auto& record : query) {
            EXPECT_FALSE(record.IsSegment());
            record.header_ = segmentHeader;
            record.ReadGroup(readGroups[0]);
            EXPECT_TRUE(record.IsSegment());
            segmentWriter.Write(record);
        }
    }
    {
        // verify round trip was correct
        BAM::BamFile file{tempSegmentBamFn};
        const BAM::BamHeader& originalHeader = file.Header();
        auto readGroups = originalHeader.ReadGroups();
        ASSERT_EQ(readGroups.size(), 1);
        EXPECT_EQ(readGroups[0].Id(), "e51ee4ef");
        EXPECT_TRUE(readGroups[0].IsSegment());

        BAM::EntireFileQuery query{SegmentReadTests::BasicSegmentBamFn};
        for (auto& record : query) {
            EXPECT_TRUE(record.IsSegment());
        }
    }

    // clean up temp files
    std::ignore = std::remove(tempCcsBamFn.c_str());
    std::ignore = std::remove(tempSegmentBamFn.c_str());
}
