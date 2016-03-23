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

// Author: Derek Barnett, Lance Hepler

#ifdef PBBAM_TESTING
#define private public
#endif

#include <gtest/gtest.h>
#include <pbbam/ReadGroupInfo.h>
#include <vector>
using namespace PacBio::BAM;
using namespace std;

TEST(ReadGroupInfoTest, IdFromMovieNameAndReadType)
{
    ReadGroupInfo rg("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0", "HQREGION");
    EXPECT_EQ("00082ba1", rg.Id());
}

TEST(ReadGroupInfoTest, FrameCodecSetOk)
{
    ReadGroupInfo rg("test");
    rg.IpdCodec(FrameCodec::V1);
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::IPD));
    EXPECT_EQ("ip", rg.BaseFeatureTag(BaseFeature::IPD));
    EXPECT_EQ(FrameCodec::V1, rg.IpdCodec());
}

TEST(ReadGroupInfoTest, SequencingChemistryOk)
{
    using std::string;
    using std::vector;

    { // P6-C4
        const vector<string> bindingKits { "100356300", "100372700" };
        const vector<string> versions { "2.1", "2.3" };
        ReadGroupInfo rg("P6C4");
        rg.SequencingKit("100356200");
        for (const string& bk : bindingKits) {
            rg.BindingKit(bk);
            for (const string& ver : versions) {
                rg.BasecallerVersion(ver);
                EXPECT_EQ("P6-C4", rg.SequencingChemistry());
            }
        }
    }

    { // S/P1-C1/beta
        const vector<string> sequencingKits { "100-619-400", "100-711-600" };
        ReadGroupInfo rg("SP1C1");
        rg.BindingKit("100-619-300");
        rg.BasecallerVersion("3.0");
        for (const string& sk : sequencingKits) {
            rg.SequencingKit(sk);
            EXPECT_EQ("S/P1-C1/beta", rg.SequencingChemistry());
        }
    }

    // basecaller 3.1.x
    { 
        const vector<string> sequencingKits { "100-619-400", "100-711-600", "100-620-000" };
        ReadGroupInfo rg("3.1");
        rg.BindingKit("100-619-300");
        rg.BasecallerVersion("3.1.0.171835");
        for (const string& sk : sequencingKits) {
            rg.SequencingKit(sk);
            EXPECT_EQ("S/P1-C1", rg.SequencingChemistry());
        }
    }
}

TEST(ReadGroupInfoTest, SequencingChemistryThrowsOnBadTriple)
{
    try {
        ReadGroupInfo rg("BAD");
        rg.BindingKit("100372700");
        rg.SequencingKit("100-619-400");
        rg.BasecallerVersion("2.0");
        //EXPECT_THROW(rg.SequencingChemistry(), InvalidSequencingChemistryException);
    } catch (InvalidSequencingChemistryException& e) {
        EXPECT_EQ(string("100372700"),   e.BindingKit());
        EXPECT_EQ(string("100-619-400"), e.SequencingKit());
        EXPECT_EQ(string("2.0"),         e.BasecallerVersion());
    }
}




