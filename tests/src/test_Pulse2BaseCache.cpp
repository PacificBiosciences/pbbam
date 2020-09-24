// Author: Derek Barnett

#include <pbbam/../../src/Pulse2BaseCache.h>

#include <gtest/gtest.h>

TEST(BAM_Pulse2BaseCache, can_determine_pulse_counts)
{
    const std::string pulseCalls{"ACccTTAGtTCAtG"};
    const std::string trimmedPC{"ACTTAGTCAG"};

    const PacBio::BAM::Pulse2BaseCache cache{pulseCalls};
    EXPECT_EQ(pulseCalls.size(), cache.NumPulses());
    EXPECT_EQ(trimmedPC.size(), cache.NumBases());
}

TEST(BAM_Pulse2BaseCache, can_remove_squashed_pulses_from_string)
{
    const std::string pulseCalls{"ACccTTAGtTCAtG"};
    const std::string trimmedPC{"ACTTAGTCAG"};
    const std::string altLabel{"-G--A--T--AC--"};
    const std::string trimmedAlt{"-GA--T-AC-"};

    const PacBio::BAM::Pulse2BaseCache cache{pulseCalls};
    EXPECT_EQ(trimmedPC, cache.RemoveSquashedPulses(pulseCalls));
    EXPECT_EQ(trimmedAlt, cache.RemoveSquashedPulses(altLabel));
}

TEST(BAM_Pulse2BaseCache, can_remove_squashed_pulses_from_integer_vector)
{
    const std::string pulseCalls{"ACccTTAGtTCAtG"};
    const std::vector<uint16_t> pkMean{5, 4, 2, 2, 3, 8, 8, 8, 4, 7, 7, 7, 3, 4};
    const std::vector<uint16_t> trimmedPkmean{5, 4, 3, 8, 8, 8, 7, 7, 7, 4};

    const PacBio::BAM::Pulse2BaseCache cache{pulseCalls};
    EXPECT_EQ(trimmedPkmean, cache.RemoveSquashedPulses(pkMean));
}
