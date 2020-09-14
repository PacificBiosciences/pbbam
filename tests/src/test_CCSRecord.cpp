// Author: Derek Barnett

#include <pbbam/ccs/CCSRecord.h>

#include <string>

#include <gtest/gtest.h>

using namespace PacBio;

TEST(CCS_CCSRecord, can_convert_to_read)
{
    const int32_t holeNumber = 77;
    const Data::Position qStart = 1000;
    const Data::Position qEnd = 1010;
    const Data::LocalContextFlags ctxtFlags =
        Data::LocalContextFlags::ADAPTER_BEFORE | Data::LocalContextFlags::ADAPTER_AFTER;
    const Data::Accuracy acc = 0.95f;
    const Data::SNR snr{0.4, 0.4, 0.4, 0.4};
    const std::string seq{"GGTTAACCAA"};
    const Data::Frames pw{3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    const std::string movie = "movie";
    const std::string chemistry = "chemistry";

    const CCS::CCSRecord ccsRecord{holeNumber, qStart, qEnd, ctxtFlags, acc, snr, seq, pw};

    const auto read = ccsRecord.ToRead(movie, chemistry);
    EXPECT_EQ(holeNumber, read.Id.HoleNumber);
    EXPECT_EQ(qStart, read.QueryStart);
    EXPECT_EQ(qEnd, read.QueryEnd);
    EXPECT_EQ(ctxtFlags, read.Flags);
    EXPECT_EQ(acc, read.ReadAccuracy);
    EXPECT_EQ(snr, read.SignalToNoise);
    EXPECT_EQ(seq, read.Seq);
    EXPECT_EQ(pw, read.PulseWidth);
}
