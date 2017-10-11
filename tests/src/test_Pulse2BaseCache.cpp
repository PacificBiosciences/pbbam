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

#include <cstdint>
#include <string>

#include <gtest/gtest.h>

#define private public

#include "PbbamTestData.h"

#include <pbbam/../../src/Pulse2BaseCache.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

TEST(Pulse2BaseCacheTest, CountsDetectedInConstructor)
{
    const string pulseCalls = "ACccTTAGtTCAtG";
    const string trimmedPC  = "ACTTAGTCAG";

    const Pulse2BaseCache cache{ pulseCalls };

    EXPECT_EQ(pulseCalls.size(), cache.NumPulses());
    EXPECT_EQ(trimmedPC.size(),  cache.NumBases());
}

TEST(Pulse2BaseCacheTest, RemovesSquashedPulsesFromString)
{
    const string pulseCalls = "ACccTTAGtTCAtG";
    const string trimmedPC  = "ACTTAGTCAG";
    const string altLabel   = "-G--A--T--AC--";
    const string trimmedAlt = "-GA--T-AC-";

    const Pulse2BaseCache cache{ pulseCalls };

    EXPECT_EQ(trimmedPC,  cache.RemoveSquashedPulses(pulseCalls));
    EXPECT_EQ(trimmedAlt, cache.RemoveSquashedPulses(altLabel));
}

TEST(Pulse2BaseCacheTest, RemovesSquashedPulsesFromVector)
{
    const string pulseCalls = "ACccTTAGtTCAtG";
    const vector<uint16_t> pkMean        = {5,4,2,2,3,8,8,8,4,7,7,7,3,4};
    const vector<uint16_t> trimmedPkmean = {5,4,3,8,8,8,7,7,7,4};

    const Pulse2BaseCache cache{ pulseCalls };

    EXPECT_EQ(trimmedPkmean, cache.RemoveSquashedPulses(pkMean));
}
