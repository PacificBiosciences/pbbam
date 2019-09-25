// Author: Derek Barnett, Lance Hepler

#include <cstddef>
#include <cstdlib>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/ReadGroupInfo.h>
#include <pbbam/exception/BundleChemistryMappingException.h>
#include <pbbam/exception/InvalidSequencingChemistryException.h>

// clang-format off

using namespace PacBio::BAM;

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
    {   // S/P3-C3/5.0 (Release 6.0)
        const std::string chem{"S/P3-C3/5.0"};
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-500-400", "101-427-500", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-500-400", "101-427-800", "5.0"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("101-500-400")
          .SequencingKit("101-427-500")
          .BasecallerVersion("5.0");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }

    {   // S/P3-C1/5.0-8M (Release 7.0)
        const std::string chem{"S/P3-C1/5.0-8M"};
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-490-800", "101-644-500", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-490-800", "101-717-100", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-717-300", "101-644-500", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-717-300", "101-717-100", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-717-400", "101-644-500", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-717-400", "101-717-100", "5.0"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("101-490-800")
          .SequencingKit("101-644-500")
          .BasecallerVersion("5.0");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }

    {   // S/P4-C2/5.0-8M (Release 8.0)
        const std::string chem{"S/P4-C2/5.0-8M"};
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-789-500", "101-826-100", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-789-500", "101-820-300", "5.0"));

        ReadGroupInfo rg("dummy");
        rg.BindingKit("101-789-500")
          .SequencingKit("101-826-100")
          .BasecallerVersion("5.0");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }
    {   // S/P1-C1.3 (Goat)
        const std::string chem{"S/P1-C1.3"};
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
    EXPECT_EQ(0, setenv(varname, PbbamTestsConfig::Data_Dir.c_str(), 0));

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
        ReadGroupInfo rg2("BAD");
        rg2.BindingKit("100372700")
          .SequencingKit("100-619-400")
          .BasecallerVersion("2.0");
    } catch (InvalidSequencingChemistryException& e) {
        EXPECT_EQ(std::string("100372700"),   e.BindingKit());
        EXPECT_EQ(std::string("100-619-400"), e.SequencingKit());
        EXPECT_EQ(std::string("2.0"),         e.BasecallerVersion());
    }
}

TEST(ReadGroupInfoTest, BasecallerVersion)
{
    // too short
    {
        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-867-300")
          .BasecallerVersion("3");
        EXPECT_THROW(rg.SequencingChemistry(), std::exception);
    }

    // initial implementation assumed single digit version numbers:
    //    const std::string ver{ basecallerVersion.substr(0, 3) };
    // So '3.299.dummy' would incorrectly be interpreted as (OK) '3.2'.
    // 3.

    try {
        ReadGroupInfo rg("dummy");
        rg.BindingKit("100-619-300")
          .SequencingKit("100-867-300")
          .BasecallerVersion("3.199.dummy");
        const std::string chem = rg.SequencingChemistry();
//        ()chem;

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
    rg.BaseFeatureTag(BaseFeature::PULSE_EXCLUSION, "pe");

    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::DELETION_QV));
    EXPECT_EQ("dq", rg.BaseFeatureTag(BaseFeature::DELETION_QV));

    rg.RemoveBaseFeature(BaseFeature::DELETION_QV);
    EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::DELETION_QV));

    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::DELETION_TAG));
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::INSERTION_QV));
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::MERGE_QV));
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::SUBSTITUTION_QV));
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::PULSE_EXCLUSION));
}

TEST(ReadGroupInfoTest, BaseIdFromBarcodedId)
{
    const ReadGroupInfo rg{"00082ba1/0--1"};
    EXPECT_EQ("00082ba1/0--1", rg.Id());
    EXPECT_EQ("00082ba1", rg.BaseId());
}

TEST(ReadGroupInfoTest, BaseIdFromNonBarcodedId)
{
    const ReadGroupInfo rg{"00082ba1"};
    EXPECT_EQ("00082ba1", rg.Id());
    EXPECT_EQ("00082ba1", rg.BaseId());
}

TEST(ReadGroupInfoTest, BarcodeDataFromBarcodedId)
{
    const ReadGroupInfo rg{"00082ba1/0--1"};
    EXPECT_EQ("00082ba1/0--1", rg.Id());
    EXPECT_EQ("00082ba1", rg.BaseId());

    const auto barcodes = rg.Barcodes();
    ASSERT_TRUE(barcodes);
    EXPECT_EQ(0, barcodes->first);
    EXPECT_EQ(1, barcodes->second);
    EXPECT_EQ(0, rg.BarcodeForward().get());
    EXPECT_EQ(1, rg.BarcodeReverse().get());
}

TEST(ReadGroupInfoTest, BarcodeDataFromIdPlusBarcodesCtor)
{
    const ReadGroupInfo rg{"00082ba1", std::pair<uint16_t, uint16_t>(0,1)};

    EXPECT_EQ("00082ba1/0--1", rg.Id());
    EXPECT_EQ("00082ba1", rg.BaseId());

    const auto barcodes = rg.Barcodes();
    ASSERT_TRUE(barcodes);
    EXPECT_EQ(0, barcodes->first);
    EXPECT_EQ(1, barcodes->second);
    EXPECT_EQ(0, rg.BarcodeForward().get());
    EXPECT_EQ(1, rg.BarcodeReverse().get());
}

TEST(ReadGroupInfoTest, NoBarcodeDataFromNonbarcodedId)
{
    {   // "standard" ID
        const ReadGroupInfo rg{"00082ba1"};
        EXPECT_EQ("00082ba1", rg.Id());
        EXPECT_EQ("00082ba1", rg.BaseId());

        const auto barcodes = rg.Barcodes();
        EXPECT_EQ(boost::none, barcodes);
        EXPECT_EQ(boost::none, rg.BarcodeForward());
        EXPECT_EQ(boost::none, rg.BarcodeReverse());
    }
    {   // no '/' found
        const ReadGroupInfo rg{"00082ba1.0--1"};
        const auto barcodes = rg.Barcodes();
        EXPECT_EQ(boost::none, barcodes);
        EXPECT_EQ(boost::none, rg.BarcodeForward());
        EXPECT_EQ(boost::none, rg.BarcodeReverse());
    }
}

TEST(ReadGroupInfoTest, NoBarcodeDataFromEmptyId)
{
    const ReadGroupInfo rg{""};
    const auto barcodes = rg.Barcodes();
    EXPECT_EQ(boost::none, barcodes);
    EXPECT_EQ(boost::none, rg.BarcodeForward());
    EXPECT_EQ(boost::none, rg.BarcodeReverse());
}

TEST(ReadGroupInfoTest, ThrowsOnConstructingIdFromMalformattedBarcodeLabels)
{
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0-1"}, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0---1"}, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0..1"};, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0"}, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/A--B"}, std::runtime_error);
}

// clang-format on
