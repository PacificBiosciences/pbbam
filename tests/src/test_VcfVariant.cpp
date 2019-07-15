// Author: Derek Barnett

#include <gtest/gtest.h>
#include <pbbam/vcf/VcfVariant.h>

using InfoField = PacBio::VCF::InfoField;
using VcfVariant = PacBio::VCF::VcfVariant;

namespace VcfVariantTests {

static const std::string BasicVariantText{
    "chrXVI\t660831\tpbsv.INS.21\tC\tCAAAGGAATGGTAAAGATGGGGGGTCAACGGACAAGGGAAAGGATCCATGGGGGCA\t."
    "\tPASS"
    "\tIMPRECISE;SVTYPE=INS;END=660831;SVLEN=55;MULTI=1,2,3\tGT:AD:DP:AC\t0/1:2:5:1,2"};

}  // namespace VcfVariantTests

TEST(VCF_Variant, default_ctor_provides_proper_default_values)
{
    VcfVariant v;

    EXPECT_TRUE(v.Chrom().empty());
    EXPECT_EQ(PacBio::BAM::UnmappedPosition, v.Position());
    EXPECT_TRUE(v.Id().empty());
    EXPECT_TRUE(v.RefAllele().empty());
    EXPECT_TRUE(v.AltAllele().empty());
    EXPECT_TRUE(v.IsQualityMissing());
    EXPECT_EQ("PASS", v.Filter());

    EXPECT_FALSE(v.IsDeletion());
    EXPECT_FALSE(v.IsInsertion());
    EXPECT_FALSE(v.IsSnp());
}

TEST(VCF_Variant, can_create_snp)
{
    const VcfVariant v{"var_snp", "3", 3000, "C", "G"};

    EXPECT_EQ("3", v.Chrom());
    EXPECT_EQ(3000, v.Position());
    EXPECT_EQ("var_snp", v.Id());
    EXPECT_EQ("C", v.RefAllele());
    EXPECT_EQ("G", v.AltAllele());
    EXPECT_TRUE(v.IsQualityMissing());
    EXPECT_EQ("PASS", v.Filter());

    EXPECT_FALSE(v.IsDeletion());
    EXPECT_FALSE(v.IsInsertion());
    EXPECT_TRUE(v.IsSnp());
}

TEST(VCF_Variant, can_create_insertion)
{
    const VcfVariant v{"var_ins", "3", 3000, "C", "CTAG"};

    EXPECT_EQ("3", v.Chrom());
    EXPECT_EQ(3000, v.Position());
    EXPECT_EQ("var_ins", v.Id());
    EXPECT_EQ("C", v.RefAllele());
    EXPECT_EQ("CTAG", v.AltAllele());
    EXPECT_TRUE(v.IsQualityMissing());
    EXPECT_EQ("PASS", v.Filter());

    EXPECT_FALSE(v.IsDeletion());
    EXPECT_TRUE(v.IsInsertion());
    EXPECT_FALSE(v.IsSnp());
}

TEST(VCF_Variant, can_create_deletion)
{
    const VcfVariant v{"var_del", "3", 3000, "TCG", "T"};

    EXPECT_EQ("3", v.Chrom());
    EXPECT_EQ(3000, v.Position());
    EXPECT_EQ("var_del", v.Id());
    EXPECT_EQ("TCG", v.RefAllele());
    EXPECT_EQ("T", v.AltAllele());
    EXPECT_TRUE(v.IsQualityMissing());
    EXPECT_EQ("PASS", v.Filter());

    EXPECT_TRUE(v.IsDeletion());
    EXPECT_FALSE(v.IsInsertion());
    EXPECT_FALSE(v.IsSnp());
}

TEST(VCF_Variant, can_determine_if_info_field_is_present)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    EXPECT_TRUE(v.HasInfoField("SVLEN"));
    EXPECT_FALSE(v.HasInfoField("nope"));
}

TEST(VCF_Variant, can_fetch_single_value_info_field)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    const auto& value = v.InfoValue("SVTYPE");
    EXPECT_TRUE(value.is_initialized());
    EXPECT_EQ("INS", value.get());
}

TEST(VCF_Variant, can_add_single_value_info_field)
{
    VcfVariant v{VcfVariantTests::BasicVariantText};

    InfoField i;
    i.id = "NEW";
    i.value = "42";
    v.AddInfoField(i);

    EXPECT_TRUE(v.HasInfoField("NEW"));
    EXPECT_EQ("42", v.InfoValue("NEW").get());
}

TEST(VCF_Variant, can_fetch_multi_value_info_field)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    const auto& values = v.InfoValues("MULTI");
    EXPECT_TRUE(values.is_initialized());
    EXPECT_EQ(3, values->size());
    EXPECT_EQ("1", values->at(0));
    EXPECT_EQ("2", values->at(1));
    EXPECT_EQ("3", values->at(2));
}

TEST(VCF_Variant, can_edit_single_value_info_field)
{
    VcfVariant v{VcfVariantTests::BasicVariantText};
    auto value = v.InfoValue("SVTYPE");
    EXPECT_TRUE(value.is_initialized());
    EXPECT_EQ("INS", value.get());

    v.InfoValue("SVTYPE", std::string{"FOO"});

    value = v.InfoValue("SVTYPE");
    EXPECT_TRUE(value.is_initialized());
    EXPECT_EQ("FOO", value.get());
}

TEST(VCF_Variant, can_edit_multi_value_info_field)
{
    VcfVariant v{VcfVariantTests::BasicVariantText};

    auto values = v.InfoValues("MULTI");
    EXPECT_TRUE(values.is_initialized());
    EXPECT_EQ(3, values->size());
    EXPECT_EQ("1", values->at(0));
    EXPECT_EQ("2", values->at(1));
    EXPECT_EQ("3", values->at(2));

    std::vector<std::string> newData{"42", "42", "42"};
    v.InfoValues("MULTI", newData);

    values = v.InfoValues("MULTI");
    EXPECT_TRUE(values.is_initialized());
    EXPECT_EQ(3, values->size());
    EXPECT_EQ("42", values->at(0));
    EXPECT_EQ("42", values->at(1));
    EXPECT_EQ("42", values->at(2));
}

TEST(VCF_Variant, can_add_multi_value_info_field)
{
    VcfVariant v{VcfVariantTests::BasicVariantText};
    InfoField i;
    i.id = "NEW";
    i.values = std::vector<std::string>{"42", "42", "42"};
    v.AddInfoField(i);

    EXPECT_TRUE(v.HasInfoField("NEW"));
    const auto& values = v.InfoValues("NEW");
    EXPECT_EQ(3, values->size());
    EXPECT_EQ("42", values->at(0));
    EXPECT_EQ("42", values->at(1));
    EXPECT_EQ("42", values->at(2));
}

TEST(VCF_Variant, can_remove_info_field)
{
    VcfVariant v{VcfVariantTests::BasicVariantText};

    EXPECT_TRUE(v.HasInfoField("SVLEN"));
    EXPECT_EQ("INS", v.InfoValue("SVTYPE").get());

    v.RemoveInfoField("SVLEN");

    EXPECT_FALSE(v.HasInfoField("SVLEN"));
    EXPECT_EQ("INS", v.InfoValue("SVTYPE").get());
}

TEST(VCF_Variant, can_fetch_all_genotype_ids)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    const auto& genotypeIds = v.GenotypeIds();
    ASSERT_EQ(4, genotypeIds.size());
    EXPECT_EQ("GT", genotypeIds.at(0));
    EXPECT_EQ("AD", genotypeIds.at(1));
    EXPECT_EQ("DP", genotypeIds.at(2));
    EXPECT_EQ("AC", genotypeIds.at(3));
}

TEST(VCF_Variant, can_fetch_all_genotype_fields)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    const auto& genotypeFields = v.Genotypes();
    ASSERT_EQ(1, genotypeFields.size());
}

TEST(VCF_Variant, can_fetch_single_value_genotype_field)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    const auto& value = v.GenotypeValue(0, "AD");
    EXPECT_TRUE(value.is_initialized());
    EXPECT_EQ("2", value.get());
}

TEST(VCF_Variant, can_fetch_multi_value_genotype_field)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    const auto& values = v.GenotypeValues(0, "AC");
    EXPECT_TRUE(values.is_initialized());
    ASSERT_EQ(2, values->size());
}

TEST(VCF_Variant, can_determine_if_sample_is_heterozygous)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    EXPECT_TRUE(v.IsSampleHeterozygous(0));
}

TEST(VCF_Variant, can_determine_if_sample_is_phased)
{
    const VcfVariant v{VcfVariantTests::BasicVariantText};
    EXPECT_FALSE(v.IsSamplePhased(0));
}
