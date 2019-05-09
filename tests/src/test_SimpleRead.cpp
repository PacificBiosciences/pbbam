// Author: Derek Barnett

#include <algorithm>

#include <gtest/gtest.h>

#include <pbbam/SimpleRead.h>
#include "../../src/Clipping.h"
#include "../../src/SimpleReadImpl.h"

using namespace PacBio::BAM;

TEST(SimpleReadTest, ClipSimpleRead)
{
    const std::string seq = "AACCGTTAGC";
    const QualityValues quals = QualityValues::FromFastq("0123456789");
    const SNR snr{0.9, 0.9, 0.9, 0.9};
    const Position qStart = 500;
    const Position qEnd = 510;
    const std::vector<uint16_t> pw{10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

    // ClipToQuery(read, 502, 509);
    const size_t clipStart = 502;
    const size_t clipEnd = 509;
    const internal::ClipResult clipResult{2, 502, 509};

    SimpleRead read{"name", seq, quals, snr, qStart, qEnd, pw};
    internal::ClipSimpleRead(read, clipResult, clipStart, clipEnd);

    const std::string expectedSeq = "CCGTTAG";
    const QualityValues expectedQuals = QualityValues::FromFastq("2345678");
    const SNR expectedSnr{0.9, 0.9, 0.9, 0.9};
    const Position expectedQStart = 502;
    const Position expectedQEnd = 509;
    const std::vector<uint16_t> expectedPw{30, 40, 50, 60, 70, 80, 90};

    EXPECT_EQ("name", read.Name);
    EXPECT_EQ(expectedSeq, read.Sequence);
    EXPECT_EQ(expectedQuals, read.Qualities);
    EXPECT_EQ(expectedQStart, read.QueryStart);
    EXPECT_EQ(expectedQEnd, read.QueryEnd);
    EXPECT_EQ(expectedSnr, read.SignalToNoise);
    EXPECT_TRUE(read.PulseWidths);
    EXPECT_TRUE(std::equal(expectedPw.cbegin(), expectedPw.cend(), read.PulseWidths->begin()));
}

TEST(SimpleReadTest, ClipMappedSimpleRead)
{
    const std::string seq = "AACCGTTAGC";
    const QualityValues quals = QualityValues::FromFastq("0123456789");
    const SNR snr{0.9, 0.9, 0.9, 0.9};
    const Position qStart = 500;
    const Position qEnd = 510;
    const std::vector<uint16_t> pw{10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    const Strand strand = Strand::FORWARD;
    const Position tStart = 100;
    const Position tEnd = 111;
    const Cigar cigar{"4=1D2I2D4="};
    const uint8_t mapQV = 80;

    // ClipToReference(read, 102, 107);
    const size_t clipStart = 102;
    const size_t clipEnd = 107;
    const internal::ClipResult clipResult{2, 502, 507, 102, Cigar{"2=1D2I2D"}};

    MappedSimpleRead read{
        SimpleRead{"name", seq, quals, snr, qStart, qEnd, pw}, strand, tStart, tEnd, cigar, mapQV};
    internal::ClipMappedRead(read, clipResult, clipStart, clipEnd);

    const std::string expectedSeq = "CCGTT";
    const QualityValues expectedQuals = QualityValues::FromFastq("23456");
    const SNR expectedSnr{0.9, 0.9, 0.9, 0.9};
    const Position expectedQStart = 502;
    const Position expectedQEnd = 507;
    const std::vector<uint16_t> expectedPw{30, 40, 50, 60, 70};
    const Strand expectedStrand = Strand::FORWARD;
    const Position expectedTStart = 102;
    const Position expectedTEnd = 107;
    const Cigar expectedCigar{"2=1D2I2D"};
    const uint8_t expectedMapQV = 80;

    EXPECT_EQ("name", read.Name);
    EXPECT_EQ(expectedSeq, read.Sequence);
    EXPECT_EQ(expectedQuals, read.Qualities);
    EXPECT_EQ(expectedQStart, read.QueryStart);
    EXPECT_EQ(expectedQEnd, read.QueryEnd);
    EXPECT_EQ(expectedSnr, read.SignalToNoise);
    EXPECT_TRUE(std::equal(expectedPw.cbegin(), expectedPw.cend(), read.PulseWidths->begin()));
    EXPECT_EQ(expectedStrand, read.Strand);
    EXPECT_EQ(expectedTStart, read.TemplateStart);
    EXPECT_EQ(expectedTEnd, read.TemplateEnd);
    EXPECT_EQ(expectedCigar, read.Cigar);
    EXPECT_EQ(expectedMapQV, read.MapQuality);
}
