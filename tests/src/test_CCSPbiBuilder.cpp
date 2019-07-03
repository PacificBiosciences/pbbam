// Author: Derek Barnett

#include <cstdio>

#include <gtest/gtest.h>

#include <pbbam/PbiRawData.h>
#include <pbbam/ccs/CCSPbiBuilder.h>

using namespace PacBio;
using namespace PacBio::CCS;

using Frames = PacBio::BAM::Frames;
using LocalContextFlags = PacBio::BAM::LocalContextFlags;

// clang-format off

namespace CCSPbiBuilderTests {

const CCSRecord& ValidRecord()
{
    static const CCSRecord record = [](){
        CCSRecord r;
        r.HoleNumber = 4391137;
        r.QueryStart = 0;
        r.QueryEnd = 459;
        r.LocalContextFlags = LocalContextFlags::ADAPTER_AFTER;
        r.Accuracy = 0.8f;
        r.SignalToNoise = {7.6, 13.9, 7.0, 12.2};
        r.Sequence = "GATTACA";
        r.PulseWidths = Frames{std::vector<uint16_t>{13, 8, 3, 14, 18, 3}};
        return r;
    }();
    return record;
}

}  // namespace CCSPbiBuilderTests

TEST(CCSPbiBuilderTest, can_create_pbi_file_from_ccs_records)
{
    const std::string pbiFilename{"test.pbi"};

    {
        const auto& record = CCSPbiBuilderTests::ValidRecord();
        CCSPbiBuilder builder{pbiFilename, "test"};
        EXPECT_EQ("test", builder.MovieName());
        builder.AddRecord(record);
        builder.AddRecord(record);
        builder.AddRecord(record);
        builder.Close();
    }
    {
        PacBio::BAM::PbiRawData index{pbiFilename};
        ASSERT_EQ(3, index.NumReads());

        const auto& basicData = index.BasicData();
        EXPECT_EQ(1610789639, basicData.rgId_[0]);
        EXPECT_EQ(1610789639, basicData.rgId_[1]);
        EXPECT_EQ(1610789639, basicData.rgId_[2]);

        EXPECT_EQ(4391137, basicData.holeNumber_[0]);
        EXPECT_EQ(4391137, basicData.holeNumber_[1]);
        EXPECT_EQ(4391137, basicData.holeNumber_[2]);

        EXPECT_EQ(0, basicData.qStart_[0]);
        EXPECT_EQ(0, basicData.qStart_[1]);
        EXPECT_EQ(0, basicData.qStart_[2]);

        EXPECT_EQ(459, basicData.qEnd_[0]);
        EXPECT_EQ(459, basicData.qEnd_[1]);
        EXPECT_EQ(459, basicData.qEnd_[2]);
    }

    remove(pbiFilename.c_str());
}

// clang-format on