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

#include "TestData.h"
#include <gtest/gtest.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/exception/BundleChemistryMappingException.h>
#include <pbbam/exception/InvalidSequencingChemistryException.h>
#include <string>
#include <vector>
#include <cstdlib>
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
    {   // P6-C4
        const string& chem = "P6-C4";
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100356300","100356200","2.1"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100356300","100356200","2.3"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100356300","100612400","2.1"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100356300","100612400","2.3"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100372700","100356200","2.1"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100372700","100356200","2.3"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100372700","100612400","2.1"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100372700","100612400","2.3"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("100356300")
          .SequencingKit("100356200")
          .BasecallerVersion("2.1");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }

    {   // S/P1-C1/beta
        const string& chem = "S/P1-C1/beta";
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-620-000","3.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-620-000","3.1"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-620-000")
          .BasecallerVersion("3.0");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }

    {   // S/P1-C1.1 (Echidna)
        const string& chem = "S/P1-C1.1";
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-867-300","3.1"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-867-300","3.2"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-867-300","3.3"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-867-300")
          .BasecallerVersion("3.1");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }

    {   // S/P1-C1.2 (Flea)
        const string& chem = "S/P1-C1.2";
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-902-100","3.1"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-902-100","3.2"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-902-100","3.3"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-902-100")
          .BasecallerVersion("3.1");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }
    {   // S/P1-C1.3 (Goat)
        const string& chem = "S/P1-C1.3";
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-972-200","3.2"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("100-619-300","100-972-200","3.3"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-972-200")
          .BasecallerVersion("3.3");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }
}

#ifdef _WIN32
int setenv(const char* name, const char* value, int overwrite)
{
    int err = 0;
    if (!overwrite) {
        size_t sz = 0;
        err = getenv_s(&sz, NULL, 0, name);
        if (err || sz) return err;
    }
    return _putenv_s(name, value);
}

int unsetenv(const char* name) {
    static const char* empty = "";
    return _putenv_s(name, empty);
}
#endif

TEST(ReadGroupInfoTest, SequencingChemistryFromMappingXml)
{
    ReadGroupInfo rg("MAYBE");
    rg.BindingKit("1").SequencingKit("2").BasecallerVersion("3.4");
    EXPECT_THROW(rg.SequencingChemistry(), InvalidSequencingChemistryException);

    // set the magic environment variable
    const char* varname = "SMRT_CHEMISTRY_BUNDLE_DIR";
    EXPECT_EQ(0, setenv(varname, tests::Data_Dir.c_str(), 0));

    EXPECT_EQ("FOUND", rg.SequencingChemistry());

    // unset the environment variable
    EXPECT_EQ(0, unsetenv(varname));

    // test memoization
    EXPECT_THROW(ReadGroupInfo::SequencingChemistryFromTriple("1", "2", "3.4"),
                 InvalidSequencingChemistryException);
    EXPECT_EQ("FOUND", rg.SequencingChemistry());

    EXPECT_EQ(0, setenv(varname, "/dev/null", 0));

    // test that a bogus SMRT_CHEMISTRY_BUNDLE_DIR throws
    EXPECT_THROW(ReadGroupInfo::SequencingChemistryFromTriple("1", "2", "3.4"),
                 BundleChemistryMappingException);

    EXPECT_EQ(0, unsetenv(varname));
}

TEST(ReadGroupInfoTest, SequencingChemistryThrowsOnBadTriple)
{
    // check that we actually throw
    ReadGroupInfo rg("BAD");
    rg.BindingKit("100372700")
      .SequencingKit("100-619-400")
      .BasecallerVersion("2.0");
    EXPECT_THROW(rg.SequencingChemistry(), InvalidSequencingChemistryException);

    // now check thrown contents
    try {
        ReadGroupInfo rg("BAD");
        rg.BindingKit("100372700")
          .SequencingKit("100-619-400")
          .BasecallerVersion("2.0");
    } catch (InvalidSequencingChemistryException& e) {
        EXPECT_EQ(string("100372700"),   e.BindingKit());
        EXPECT_EQ(string("100-619-400"), e.SequencingKit());
        EXPECT_EQ(string("2.0"),         e.BasecallerVersion());
    }
}

TEST(ReadGroupInfoTest, BasecallerVersion)
{
    // too short
    try {
        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-867-300")
          .BasecallerVersion("3");
        const string chem = rg.SequencingChemistry();
        (void)chem;

    } catch (std::runtime_error& e) {
        EXPECT_EQ(string("basecaller version too short: 3"), string(e.what()));
    }

    // initial implementation assumed single digit version numbers:
    //    const string ver{ basecallerVersion.substr(0, 3) };
    // So '3.299.dummy' would incorrectly be interpreted as (OK) '3.2'.
    // 3.

    try {
        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-867-300")
          .BasecallerVersion("3.199.dummy");   
        const string chem = rg.SequencingChemistry();
        (void)chem;

    } catch (InvalidSequencingChemistryException& e) {
        EXPECT_EQ("100-619-300", e.BindingKit());
        EXPECT_EQ("100-867-300", e.SequencingKit());
        EXPECT_EQ("3.199.dummy", e.BasecallerVersion());
    }
    //EXPECT_THROW(rg.SequencingChemistry(), InvalidSequencingChemistryException);
}

TEST(ReadGroupInfoTest, ClearBaseFeatures)
{
    ReadGroupInfo rg("test");
    rg.BaseFeatureTag(BaseFeature::DELETION_QV,     "dq");
    rg.BaseFeatureTag(BaseFeature::DELETION_TAG,    "dt");
    rg.BaseFeatureTag(BaseFeature::INSERTION_QV,    "iq");
    rg.BaseFeatureTag(BaseFeature::MERGE_QV,        "mq");
    rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_QV, "sq");
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::DELETION_QV));
    EXPECT_EQ("dq", rg.BaseFeatureTag(BaseFeature::DELETION_QV));

    rg.ClearBaseFeatures();
    EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::DELETION_QV));
    EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::DELETION_TAG));
    EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::INSERTION_QV));
    EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::MERGE_QV));
    EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::SUBSTITUTION_QV));
}

TEST(ReadGroupInfoTest, RemoveBaseFeature)
{
    ReadGroupInfo rg("test");
    rg.BaseFeatureTag(BaseFeature::DELETION_QV,     "dq");
    rg.BaseFeatureTag(BaseFeature::DELETION_TAG,    "dt");
    rg.BaseFeatureTag(BaseFeature::INSERTION_QV,    "iq");
    rg.BaseFeatureTag(BaseFeature::MERGE_QV,        "mq");
    rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_QV, "sq");
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::DELETION_QV));
    EXPECT_EQ("dq", rg.BaseFeatureTag(BaseFeature::DELETION_QV));

    rg.RemoveBaseFeature(BaseFeature::DELETION_QV);
    EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::DELETION_QV));
    
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::DELETION_TAG));
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::INSERTION_QV));
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::MERGE_QV));
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::SUBSTITUTION_QV));
}

