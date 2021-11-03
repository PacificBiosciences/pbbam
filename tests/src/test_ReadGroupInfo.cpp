#include <pbbam/ReadGroupInfo.h>

#include <cstddef>
#include <cstdlib>

#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include <boost/optional/optional_io.hpp>

#include <pbbam/exception/BundleChemistryMappingException.h>
#include <pbbam/exception/InvalidSequencingChemistryException.h>

#include "PbbamTestData.h"

// clang-format off

using namespace PacBio::BAM;

TEST(BAM_ReadGroupInfo, can_generate_base_id_from_id_string)
{
    const std::string rg{"123456578"};
    EXPECT_EQ("123456578", ReadGroupInfo::GetBaseId(rg));
}

TEST(BAM_ReadGroupInfo, can_generate_base_id_from_id_string_with_barcodes)
{
    const std::string rg{"123456578/0--0"};
    EXPECT_EQ("123456578", ReadGroupInfo::GetBaseId(rg));
}

TEST(BAM_ReadGroupInfo, can_generate_id_from_movie_and_read_type)
{
    const ReadGroupInfo rg{"m140905_042212_sidney_c100564852550000001823085912221377_s1_X0", "HQREGION"};
    EXPECT_EQ("00082ba1", rg.Id());
}

TEST(BAM_ReadGroupInfo, can_describe_frame_codec)
{
    ReadGroupInfo rg{"test"};
    rg.IpdCodec(FrameCodec::V1);
    EXPECT_TRUE(rg.HasBaseFeature(BaseFeature::IPD));
    EXPECT_EQ("ip", rg.BaseFeatureTag(BaseFeature::IPD));
    EXPECT_EQ(FrameCodec::V1, rg.IpdCodec());
}

TEST(BAM_ReadGroupInfo, can_lookup_chemistry_from_compiled_chemistry_table)
{
    {   // S/P3-C3/5.0 (Release 6.0)
        const std::string chem{"S/P3-C3/5.0"};
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-500-400", "101-427-500", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-500-400", "101-427-800", "5.0"));

        ReadGroupInfo rg{"dummy"};
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

        ReadGroupInfo rg{"dummy"};
        rg.BindingKit("101-490-800")
          .SequencingKit("101-644-500")
          .BasecallerVersion("5.0");
        EXPECT_EQ(chem, rg.SequencingChemistry());
    }

    {   // S/P4-C2/5.0-8M (Release 8.0)
        const std::string chem{"S/P4-C2/5.0-8M"};
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-789-500", "101-826-100", "5.0"));
        EXPECT_EQ(chem, ReadGroupInfo::SequencingChemistryFromTriple("101-789-500", "101-820-300", "5.0"));

        ReadGroupInfo rg{"dummy"};
        rg.BindingKit("101-789-500")
          .SequencingKit("101-826-100")
          .BasecallerVersion("5.0");
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

TEST(BAM_ReadGroupInfo, can_lookup_chemistry_from_mapping_xml)
{
    ReadGroupInfo rg{"MAYBE"};
    rg.BindingKit("1")
      .SequencingKit("2")
      .BasecallerVersion("3.4");
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

TEST(BAM_ReadGroupInfo, throws_on_bad_chemistry_triple)
{
    try {
        ReadGroupInfo rg{"BAD"};
        rg.BindingKit("100372700")
          .SequencingKit("100-619-400")
          .BasecallerVersion("2.0");
        const auto chem = rg.SequencingChemistry();
        ASSERT_TRUE(false);
    } catch (InvalidSequencingChemistryException& e) {
        EXPECT_EQ("100372700",   e.BindingKit());
        EXPECT_EQ("100-619-400", e.SequencingKit());
        EXPECT_EQ("2.0",         e.BasecallerVersion());
    }
}

TEST(BAM_ReadGroupInfo, throws_on_invalid_basecaller_version)
{
    // too short
    {
        ReadGroupInfo rg{"dummy"};
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
        ReadGroupInfo rg{"dummy"};
        rg.BindingKit("100-619-300")
          .SequencingKit("100-867-300")
          .BasecallerVersion("3.199.dummy");
        const auto chem = rg.SequencingChemistry();
        ASSERT_TRUE(false);
    } catch (InvalidSequencingChemistryException& e) {
        EXPECT_EQ("100-619-300", e.BindingKit());
        EXPECT_EQ("100-867-300", e.SequencingKit());
        EXPECT_EQ("3.199.dummy", e.BasecallerVersion());
    }
    //EXPECT_THROW(rg.SequencingChemistry(), InvalidSequencingChemistryException);
}

TEST(BAM_ReadGroupInfo, can_clear_all_base_features)
{
    ReadGroupInfo rg{"test"};
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

TEST(BAM_ReadGroupInfo, can_remove_single_base_feature)
{
    ReadGroupInfo rg{"test"};
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

TEST(BAM_ReadGroupInfo, can_fetch_id_types_from_barcoded_id)
{
    const ReadGroupInfo rg{"00082ba1/0--1"};
    EXPECT_EQ("00082ba1/0--1", rg.Id());
    EXPECT_EQ("00082ba1", rg.BaseId());
}

TEST(BAM_ReadGroupInfo, can_fetch_id_types_from_standard_id)
{
    const ReadGroupInfo rg{"00082ba1"};
    EXPECT_EQ("00082ba1", rg.Id());
    EXPECT_EQ("00082ba1", rg.BaseId());
}

TEST(BAM_ReadGroupInfo, can_determine_barcodes_from_barcoded_id_string)
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

TEST(BAM_ReadGroupInfo, can_determine_barcodes_from_id_string_and_barcode_pair)
{
    const ReadGroupInfo rg{"00082ba1", std::make_pair(0,1)};
    EXPECT_EQ("00082ba1/0--1", rg.Id());
    EXPECT_EQ("00082ba1", rg.BaseId());

    const auto barcodes = rg.Barcodes();
    ASSERT_TRUE(barcodes);
    EXPECT_EQ(0, barcodes->first);
    EXPECT_EQ(1, barcodes->second);
    EXPECT_EQ(0, rg.BarcodeForward().get());
    EXPECT_EQ(1, rg.BarcodeReverse().get());
}

TEST(BAM_ReadGroupInfo, returns_no_barcodes_from_non_barcoded_id)
{
    {   // "standard" ID
        const ReadGroupInfo rg{"00082ba1"};
        EXPECT_EQ("00082ba1", rg.Id());
        EXPECT_EQ("00082ba1", rg.BaseId());

        const auto barcodes = rg.Barcodes();
        EXPECT_FALSE(barcodes);
        EXPECT_EQ(boost::none, rg.BarcodeForward());
        EXPECT_EQ(boost::none, rg.BarcodeReverse());
    }
    {   // no '/' found
        const ReadGroupInfo rg{"00082ba1.0--1"};
        const auto barcodes = rg.Barcodes();
        EXPECT_FALSE(barcodes);
        EXPECT_EQ(boost::none, rg.BarcodeForward());
        EXPECT_EQ(boost::none, rg.BarcodeReverse());
    }
}

TEST(BAM_ReadGroupInfo, returns_no_barcodes_from_empty_id)
{
    const ReadGroupInfo rg{""};
    const auto barcodes = rg.Barcodes();
    EXPECT_FALSE(barcodes);
    EXPECT_EQ(boost::none, rg.BarcodeForward());
    EXPECT_EQ(boost::none, rg.BarcodeReverse());
}

TEST(BAM_ReadGroupInfo, throws_on_malformatted_barcoded_ids)
{
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0-1"}, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0---1"}, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0..1"};, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/0"}, std::runtime_error);
    EXPECT_THROW(ReadGroupInfo{"00082ba1/A--B"}, std::runtime_error);
}

TEST(BAM_ReadGroupInfo, barcodes_do_not_affect_read_group_hash)
{
    const std::string id{"d9020782"};
    const std::string movieName{"dummy_movie"};
    const std::string readType{"SUBREAD"};
    const std::pair<uint16_t, uint16_t> barcodes{32,32};

    const ReadGroupInfo rg0{id};
    const ReadGroupInfo rg1{movieName, readType};
    const ReadGroupInfo rg2{movieName, readType, barcodes};
    const ReadGroupInfo rg3{id, barcodes};

    const std::string expectedId{"d9020782"};
    const std::string expectedBarcodedId{"d9020782/32--32"};
    EXPECT_EQ(expectedId, rg0.Id());
    EXPECT_EQ(expectedId, rg1.Id());
    EXPECT_EQ(expectedBarcodedId, rg2.Id());
    EXPECT_EQ(expectedBarcodedId, rg3.Id());
}

// clang-format on
