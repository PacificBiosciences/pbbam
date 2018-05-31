// Author: Derek Barnett

#include <iostream>

#include <gtest/gtest.h>
#include <pbbam/vcf/VcfHeader.h>

using ContigDefinition = PacBio::VCF::ContigDefinition;
using FilterDefinition = PacBio::VCF::FilterDefinition;
using FormatDefinition = PacBio::VCF::FormatDefinition;
using GeneralDefinition = PacBio::VCF::GeneralDefinition;
using InfoDefinition = PacBio::VCF::InfoDefinition;
using VcfHeader = PacBio::VCF::VcfHeader;

namespace VcfHeaderTests {

static const std::string BasicHeaderText{
    "##fileformat=VCFv4.2\n"
    "##fileDate=20180509\n"
    "##contig=<ID=ctg1,length=4200,assembly=foo,md5=dead123beef>\n"
    "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant "
    "described in this record\">\n"
    "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT "
    "alleles\">\n"
    "##INFO=<ID=SVANN,Number=.,Type=String,Description=\"Repeat annotation of structural "
    "variant\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Per-sample read depth of this structural "
    "variant\">\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position for this "
    "sample\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUnnamedSample\n"};

}  // namespace VcfHeaderTests

TEST(VCF_GeneralDefinition, throws_on_missing_required_fields)
{
    const std::string id{"id"};
    const std::string desc{"desc"};

    EXPECT_THROW(GeneralDefinition("", desc), std::runtime_error);
    EXPECT_THROW(GeneralDefinition(id, ""), std::runtime_error);
}

TEST(VCF_ContigDefinition, throws_on_missing_required_fields)
{
    EXPECT_THROW(ContigDefinition(""), std::runtime_error);
}

TEST(VCF_ContigDefinition, can_edit_and_query_attributes)
{
    ContigDefinition contig{"id"};

    EXPECT_TRUE(contig.Attributes().empty());

    const std::vector<std::pair<std::string, std::string>> attributes{{"assembly", "foo"},
                                                                      {"length", "42"}};
    contig.Attributes(attributes);
    ASSERT_EQ(2, contig.Attributes().size());
    EXPECT_EQ("foo", contig.Attributes().at(0).second);
    EXPECT_EQ("42", contig.Attributes().at(1).second);

    contig.AddAttribute({"md5", "dead123beef"});
    ASSERT_EQ(3, contig.Attributes().size());
    EXPECT_EQ("dead123beef", contig.Attributes().at(2).second);
}

TEST(VCF_FilterDefinition, throws_on_missing_required_fields)
{
    const std::string id{"id"};
    const std::string desc{"desc"};

    EXPECT_THROW(FilterDefinition("", desc), std::runtime_error);
    EXPECT_THROW(FilterDefinition(id, ""), std::runtime_error);
}

TEST(VCF_InfoDefinition, throws_on_missing_required_fields)
{
    const std::string id{"id"};
    const std::string num{"num"};
    const std::string type{"type"};
    const std::string desc{"desc"};

    EXPECT_THROW(InfoDefinition("", num, type, desc), std::runtime_error);
    EXPECT_THROW(InfoDefinition(id, "", type, desc), std::runtime_error);
    EXPECT_THROW(InfoDefinition(id, num, "", desc), std::runtime_error);
    EXPECT_THROW(InfoDefinition(id, num, type, ""), std::runtime_error);
}

TEST(VCF_InfoDefinition, missing_optional_fields_is_not_error)
{
    InfoDefinition info{"id", "num", "type", "description"};

    EXPECT_FALSE(info.Source().is_initialized());
    EXPECT_FALSE(info.Version().is_initialized());

    info.Source("source");
    info.Version("version");

    EXPECT_TRUE(info.Source().is_initialized());
    EXPECT_TRUE(info.Version().is_initialized());
}

TEST(VCF_Header, defaults_to_current_version)
{
    VcfHeader hdr;
    EXPECT_EQ("VCFv4.2", hdr.Version());
}

TEST(VCF_Header, can_lookup_contig_defnition_by_id)
{
    const VcfHeader hdr{VcfHeaderTests::BasicHeaderText};
    const auto& contig = hdr.ContigDefinition("ctg1");
    ASSERT_EQ(3, contig.Attributes().size());
    EXPECT_EQ("length", contig.Attributes().at(0).first);
    EXPECT_EQ("assembly", contig.Attributes().at(1).first);
    EXPECT_EQ("md5", contig.Attributes().at(2).first);
}

TEST(VCF_Header, can_lookup_format_definition_by_id)
{
    const VcfHeader hdr{VcfHeaderTests::BasicHeaderText};
    const auto& format = hdr.FormatDefinition("GT");
    EXPECT_EQ("GT", format.Id());
}

TEST(VCF_Header, can_lookup_general_definition_by_id)
{
    const VcfHeader hdr{VcfHeaderTests::BasicHeaderText};
    const auto& def = hdr.GeneralDefinition("fileformat");
    EXPECT_EQ("fileformat", def.Id());
}

TEST(VCF_Header, can_lookup_info_definition_by_id)
{
    const VcfHeader hdr{VcfHeaderTests::BasicHeaderText};
    const auto& info = hdr.InfoDefinition("IMPRECISE");
    EXPECT_EQ("IMPRECISE", info.Id());
}

TEST(VCF_Header, can_lookup_sample)
{
    const VcfHeader hdr{VcfHeaderTests::BasicHeaderText};
    const auto idx = hdr.IndexOfSample("UnnamedSample");
    const auto sample = hdr.SampleAt(idx);
    EXPECT_EQ("UnnamedSample", sample);
}

TEST(VCF_Header, add_duplicate_format_replaces_existing_definition)
{
    VcfHeader hdr{VcfHeaderTests::BasicHeaderText};
    const auto initialFormat = hdr.FormatDefinition("GT");
    EXPECT_EQ("Genotype", initialFormat.Description());

    const FormatDefinition newFormat{"GT", "num", "type", "newDescription"};
    hdr.AddFormatDefinition(newFormat);

    const auto nowFormat = hdr.FormatDefinition("GT");
    EXPECT_EQ("newDescription", nowFormat.Description());

    // rest of defs unchanged
    const auto& formatDefs = hdr.FormatDefinitions();
    ASSERT_EQ(3, formatDefs.size());
    EXPECT_EQ("AD", formatDefs.at(1).Id());
    EXPECT_EQ("DP", formatDefs.at(2).Id());
}

TEST(VCF_Header, add_duplicate_info_replaces_existing_definition)
{
    VcfHeader hdr{VcfHeaderTests::BasicHeaderText};
    const auto initialInfo = hdr.InfoDefinition("IMPRECISE");
    EXPECT_EQ("Imprecise structural variant", initialInfo.Description());

    const InfoDefinition newInfo{"IMPRECISE", "num", "type", "newInfo"};
    hdr.AddInfoDefinition(newInfo);

    const auto nowInfo = hdr.InfoDefinition("IMPRECISE");
    EXPECT_EQ("newInfo", nowInfo.Description());

    // rest of defs unchanged
    const auto& infoDefs = hdr.InfoDefinitions();
    ASSERT_EQ(5, infoDefs.size());
    EXPECT_EQ("SVTYPE", infoDefs.at(1).Id());
    EXPECT_EQ("END", infoDefs.at(2).Id());
    EXPECT_EQ("SVLEN", infoDefs.at(3).Id());
    EXPECT_EQ("SVANN", infoDefs.at(4).Id());
}
