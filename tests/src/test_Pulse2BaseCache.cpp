// Author: Derek Barnett

#include <cstdint>
#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/../../src/Pulse2BaseCache.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

TEST(Pulse2BaseCacheTest, CountsDetectedInConstructor)
{
    const string pulseCalls = "ACccTTAGtTCAtG";
    const string trimmedPC = "ACTTAGTCAG";

    const Pulse2BaseCache cache{pulseCalls};

    EXPECT_EQ(pulseCalls.size(), cache.NumPulses());
    EXPECT_EQ(trimmedPC.size(), cache.NumBases());
}

TEST(Pulse2BaseCacheTest, RemovesSquashedPulsesFromString)
{
    const string pulseCalls = "ACccTTAGtTCAtG";
    const string trimmedPC = "ACTTAGTCAG";
    const string altLabel = "-G--A--T--AC--";
    const string trimmedAlt = "-GA--T-AC-";

    const Pulse2BaseCache cache{pulseCalls};

    EXPECT_EQ(trimmedPC, cache.RemoveSquashedPulses(pulseCalls));
    EXPECT_EQ(trimmedAlt, cache.RemoveSquashedPulses(altLabel));
}

TEST(Pulse2BaseCacheTest, RemovesSquashedPulsesFromVector)
{
    const string pulseCalls = "ACccTTAGtTCAtG";
    const vector<uint16_t> pkMean = {5, 4, 2, 2, 3, 8, 8, 8, 4, 7, 7, 7, 3, 4};
    const vector<uint16_t> trimmedPkmean = {5, 4, 3, 8, 8, 8, 7, 7, 7, 4};

    const Pulse2BaseCache cache{pulseCalls};

    EXPECT_EQ(trimmedPkmean, cache.RemoveSquashedPulses(pkMean));
}
