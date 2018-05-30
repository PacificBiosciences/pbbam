// Author: Derek Barnett

#include <gtest/gtest.h>
#include <pbbam/vcf/VcfFormat.h>
#include <pbbam/vcf/VcfHeader.h>
#include <pbbam/vcf/VcfVariant.h>

#include "PbbamTestData.h"

using ContigDefinition = PacBio::VCF::ContigDefinition;
using FilterDefinition = PacBio::VCF::FilterDefinition;
using FormatDefinition = PacBio::VCF::FormatDefinition;
using GeneralDefinition = PacBio::VCF::GeneralDefinition;
using InfoDefinition = PacBio::VCF::InfoDefinition;
using Sample = PacBio::VCF::Sample;
using VcfFormat = PacBio::VCF::VcfFormat;
using VcfHeader = PacBio::VCF::VcfHeader;
using VcfVariant = PacBio::VCF::VcfVariant;

namespace VcfFormatTests {

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
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUnnamedSample"};

// does not have ##contig line(s) in file
static const std::string FileHeaderText{
    "##fileformat=VCFv4.2\n"
    "##fileDate=20180509\n"
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
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUnnamedSample"};

static const std::string BasicVariantText{
    "chrXVI\t660831\tpbsv.INS.21\tC\tCAAAGGAATGGTAAAGATGGGGGGTCAACGGACAAGGGAAAGGATCCATGGGGGCA\t."
    "\tPASS"
    "\tIMPRECISE;SVTYPE=INS;END=660831;SVLEN=55;MULTI=1,2,3\tGT:AD:DP:AC\t0/1:2:5:1,2"};

static const std::string VcfFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                               "/vcf/structural_variants.vcf"};

}  // namespace VcfFormatTests

TEST(VCF_Format, provides_current_version)
{
    const std::string version = VcfFormat::CurrentVersion();
    EXPECT_EQ("VCFv4.2", version);
}

//## ----------------------------------------------------------------- ##
//
//              HEADER FORMATTING
//
//## ----------------------------------------------------------------- ##

TEST(VCF_Format, can_format_contig_definition)
{
    const ContigDefinition def{"ctg1",
                               {{"length", "4200"}, {"assembly", "foo"}, {"md5", "dead123beef"}}};
    const auto text = VcfFormat::FormattedContigDefinition(def);
    EXPECT_EQ("##contig=<ID=ctg1,length=4200,assembly=foo,md5=dead123beef>", text);
}

TEST(VCF_Format, can_format_filter_definition)
{
    const FilterDefinition def{"FILTER1", "Filter1"};
    const auto text = VcfFormat::FormattedFilterDefinition(def);
    EXPECT_EQ("##FILTER=<ID=FILTER1,Description=\"Filter1\">", text);
}

TEST(VCF_Format, can_format_format_definition)
{
    const FormatDefinition def{"GT", "1", "String", "Genotype"};
    const auto text = VcfFormat::FormattedFormatDefinition(def);
    EXPECT_EQ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", text);
}

TEST(VCF_Format, can_format_general_header_definition)
{
    const GeneralDefinition def{"phasing", "partial"};
    const auto text = VcfFormat::FormattedGeneralDefinition(def);
    EXPECT_EQ("##phasing=partial", text);
}

TEST(VCF_Format, can_format_info_definition)
{
    const InfoDefinition def{"IMPRECISE", "0", "Flag", "Imprecise structural variant"};
    const auto text = VcfFormat::FormattedInfoDefinition(def);
    EXPECT_EQ(
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">",
        text);
}

TEST(VCF_Format, can_format_info_definition_with_optional_fields)
{
    {  // with Source
        const InfoDefinition def{"IMPRECISE", "0", "Flag", "Imprecise structural variant",
                                 "source1"};
        const auto text = VcfFormat::FormattedInfoDefinition(def);
        EXPECT_EQ(
            "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural "
            "variant\",Source=\"source1\">",
            text);
    }

    {  // with Version
        const InfoDefinition def{"IMPRECISE", "0",       "Flag", "Imprecise structural variant",
                                 "",          "version1"};
        const auto text = VcfFormat::FormattedInfoDefinition(def);
        EXPECT_EQ(
            "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural "
            "variant\",Version=\"version1\">",
            text);
    }
    {  // with Source & Version
        const InfoDefinition def{"IMPRECISE", "0",       "Flag", "Imprecise structural variant",
                                 "source1",   "version1"};
        const auto text = VcfFormat::FormattedInfoDefinition(def);
        EXPECT_EQ(
            "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural "
            "variant\",Source=\"source1\",Version=\"version1\">",
            text);
    }
}

TEST(VCF_Format, can_format_basic_header)
{
    const VcfHeader header{VcfFormatTests::BasicHeaderText};
    const auto text = VcfFormat::FormattedHeader(header);
    EXPECT_EQ(VcfFormatTests::BasicHeaderText, text);
}

TEST(VCF_Format, format_basic_header_with_only_filedate)
{
    VcfHeader header;
    header.FileDate("1770704");
    std::string text;
    EXPECT_NO_THROW(text = VcfFormat::FormattedHeader(header));
}

TEST(VCF_Format, format_basic_header_with_only_version)
{
    VcfHeader header;
    header.Version("3.14");
    std::string text;
    EXPECT_NO_THROW(text = VcfFormat::FormattedHeader(header));
}

//## ----------------------------------------------------------------- ##
//
//              HEADER PARSING
//
//## ----------------------------------------------------------------- ##

TEST(VCF_Format, can_parse_general_header_definition)
{
    const auto phasing = VcfFormat::ParsedGeneralDefinition("##phasing=partial");
    EXPECT_EQ("phasing", phasing.Id());
    EXPECT_EQ("partial", phasing.Text());
}

TEST(VCF_Format, parsing_general_header_definition_throws_on_empty_string)
{
    EXPECT_THROW(VcfFormat::ParsedGeneralDefinition(""), std::runtime_error);
}

TEST(VCF_Format, parsing_general_header_definition_throws_on_non_vcf_input)
{
    EXPECT_THROW(VcfFormat::ParsedGeneralDefinition("not_vcf_header_line"), std::runtime_error);
    EXPECT_THROW(VcfFormat::ParsedGeneralDefinition("#line=not_vcf_header_line"),
                 std::runtime_error);
    EXPECT_THROW(VcfFormat::ParsedGeneralDefinition("##line,not_vcf_header_line"),
                 std::runtime_error);
}

TEST(VCF_Format, can_parse_contig_definition_with_id_only)
{
    const auto contig = VcfFormat::ParsedContigDefinition("##contig=<ID=ctg1>");
    EXPECT_EQ("ctg1", contig.Id());
    EXPECT_TRUE(contig.Attributes().empty());
}

TEST(VCF_Format, can_parse_contig_definition_with_attributes)
{
    const auto contig =
        VcfFormat::ParsedContigDefinition("##contig=<ID=ctg1,assembly=foo,length=3>");
    EXPECT_EQ("ctg1", contig.Id());
    ASSERT_EQ(2, contig.Attributes().size());

    const auto& firstAttr = contig.Attributes().at(0);
    EXPECT_EQ("assembly", firstAttr.first);
    EXPECT_EQ("foo", firstAttr.second);

    const auto& secondAttr = contig.Attributes().at(1);
    EXPECT_EQ("length", secondAttr.first);
    EXPECT_EQ("3", secondAttr.second);
}

TEST(VCF_Format, parsing_contig_header_definition_throws_on_malformed_contig_line)
{
    // internal code already checks for "##contig=<"

    EXPECT_THROW(VcfFormat::ParsedContigDefinition("##contig=<foo"), std::runtime_error);
    EXPECT_THROW(VcfFormat::ParsedContigDefinition("##contig=<ID=,>"), std::runtime_error);
}

TEST(VCF_Format, can_parse_filter_definition)
{
    const auto filter =
        VcfFormat::ParsedFilterDefinition("##FILTER=<ID=FILTER1,Description=\"Filter1\">\n");
    EXPECT_EQ("FILTER1", filter.Id());
    EXPECT_EQ("Filter1", filter.Description());
}

TEST(VCF_Format, parsing_filter_definition_throws_on_malformed_filter_line)
{
    // internal code already checks for "##FILTER=<"

    EXPECT_THROW(VcfFormat::ParsedFilterDefinition("##FILTER=<foo"), std::runtime_error);
    EXPECT_THROW(VcfFormat::ParsedFilterDefinition("##FILTER=<ID=,>"), std::runtime_error);
}

TEST(VCF_Format, can_parse_format_definition)
{
    const auto format = VcfFormat::ParsedFormatDefinition(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    EXPECT_EQ("GT", format.Id());
    EXPECT_EQ("1", format.Number());
    EXPECT_EQ("String", format.Type());
    EXPECT_EQ("Genotype", format.Description());
}

TEST(VCF_Format, parsing_format_definition_throws_on_malformed_filter_line)
{
    // internal code already checks for "##FORMAT=<"

    EXPECT_THROW(VcfFormat::ParsedFormatDefinition("##FORMAT=<foo"), std::runtime_error);
    EXPECT_THROW(VcfFormat::ParsedFormatDefinition("##FORMAT=<ID=,>"), std::runtime_error);
}

TEST(VCF_Format, can_parse_info_definition)
{
    const auto info = VcfFormat::ParsedInfoDefinition(
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n");
    EXPECT_EQ("IMPRECISE", info.Id());
    EXPECT_EQ("0", info.Number());
    EXPECT_EQ("Flag", info.Type());
    EXPECT_EQ("Imprecise structural variant", info.Description());
    EXPECT_FALSE(info.Source().is_initialized());
    EXPECT_FALSE(info.Version().is_initialized());
}

TEST(VCF_Format, parsing_info_definition_throws_on_malformed_info_line)
{
    // internal code already checks for "##INFO=<"

    EXPECT_THROW(VcfFormat::ParsedInfoDefinition("##INFO=<foo"), std::runtime_error);
    EXPECT_THROW(VcfFormat::ParsedInfoDefinition("##INFO=<ID=,>"), std::runtime_error);
}

TEST(VCF_Format, can_create_header_from_text)
{
    const VcfHeader hdr{VcfFormatTests::BasicHeaderText};

    EXPECT_EQ("VCFv4.2", hdr.Version());
    EXPECT_EQ("20180509", hdr.FileDate());

    const auto& infos = hdr.InfoDefinitions();
    ASSERT_EQ(5, infos.size());
    EXPECT_EQ("IMPRECISE", infos.at(0).Id());
    EXPECT_EQ("SVTYPE", infos.at(1).Id());
    EXPECT_EQ("END", infos.at(2).Id());
    EXPECT_EQ("SVLEN", infos.at(3).Id());
    EXPECT_EQ("SVANN", infos.at(4).Id());

    const auto& contigs = hdr.ContigDefinitions();
    ASSERT_EQ(1, contigs.size());
    EXPECT_EQ("ctg1", contigs.at(0).Id());

    ASSERT_EQ(3, contigs.at(0).Attributes().size());
    EXPECT_EQ("length", contigs.at(0).Attributes().at(0).first);
    EXPECT_EQ("assembly", contigs.at(0).Attributes().at(1).first);
    EXPECT_EQ("md5", contigs.at(0).Attributes().at(2).first);

    const auto& filters = hdr.FilterDefinitions();
    ASSERT_EQ(0, filters.size());

    const auto& formats = hdr.FormatDefinitions();
    ASSERT_EQ(3, formats.size());
    EXPECT_EQ("GT", formats.at(0).Id());
    EXPECT_EQ("AD", formats.at(1).Id());
    EXPECT_EQ("DP", formats.at(2).Id());

    const auto& samples = hdr.Samples();
    ASSERT_EQ(1, samples.size());
    EXPECT_EQ("UnnamedSample", samples[0]);
}

TEST(VCF_Format, header_parsing_throws_on_missing_fileformat_line)
{
    const std::string missingFormat{
        "##fileDate=20180509\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUnnamedSample\n"};

    EXPECT_THROW({ VcfHeader h(missingFormat); }, std::runtime_error);
}

TEST(VCF_Format, header_parsing_throws_on_non_vcf_header_line)
{
    const std::string nonVcfLine{
        "##fileformat=VCFv4.2\n"
        " --- how did I get in here?? --- \n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUnnamedSample\n"};

    EXPECT_THROW({ VcfHeader h(nonVcfLine); }, std::runtime_error);
}

TEST(VCF_Format, can_parse_header_from_stream)
{
    std::istringstream in(VcfFormatTests::BasicHeaderText);
    const auto header = VcfFormat::HeaderFromStream(in);
    EXPECT_EQ(VcfFormatTests::BasicHeaderText, VcfFormat::FormattedHeader(header));
}

TEST(VCF_Format, can_parse_header_from_file)
{
    const std::string fn{VcfFormatTests::VcfFn};
    const auto header = VcfFormat::HeaderFromFile(fn);
    EXPECT_EQ(VcfFormatTests::FileHeaderText, VcfFormat::FormattedHeader(header));
}

//## ----------------------------------------------------------------- ##
//
//              VARIANT FORMATTING
//
//## ----------------------------------------------------------------- ##

TEST(VCF_Format, can_format_basic_variant)
{
    const VcfVariant var = VcfFormat::ParsedVariant(VcfFormatTests::BasicVariantText);
    const auto text = VcfFormat::FormattedVariant(var);
    EXPECT_EQ(VcfFormatTests::BasicVariantText, text);
}

//## ----------------------------------------------------------------- ##
//
//              VARIANT PARSING
//
//## ----------------------------------------------------------------- ##

TEST(VCF_Format, can_create_variant_from_text)
{
    const VcfVariant var = VcfFormat::ParsedVariant(VcfFormatTests::BasicVariantText);

    // CHROM POS ID REF ALT REF QUAL FILTER
    EXPECT_EQ("chrXVI", var.Chrom());
    EXPECT_EQ(660831, var.Position());
    EXPECT_EQ("pbsv.INS.21", var.Id());
    EXPECT_EQ("C", var.RefAllele());
    EXPECT_EQ("CAAAGGAATGGTAAAGATGGGGGGTCAACGGACAAGGGAAAGGATCCATGGGGGCA", var.AltAllele());
    EXPECT_TRUE(var.IsQualityMissing());
    EXPECT_EQ("PASS", var.Filter());

    // INFO
    const auto& infoFields = var.InfoFields();
    ASSERT_EQ(5, infoFields.size());
    EXPECT_EQ("IMPRECISE", infoFields.at(0).id);
    EXPECT_EQ("SVTYPE", infoFields.at(1).id);
    EXPECT_EQ("END", infoFields.at(2).id);
    EXPECT_EQ("SVLEN", infoFields.at(3).id);
    EXPECT_EQ("MULTI", infoFields.at(4).id);

    // GENOTYPES
    const auto& ids = var.GenotypeIds();
    ASSERT_EQ(4, ids.size());
    EXPECT_EQ("GT", ids.at(0));
    EXPECT_EQ("AD", ids.at(1));
    EXPECT_EQ("DP", ids.at(2));
    EXPECT_EQ("AC", ids.at(3));

    const auto& genotypes = var.Genotypes();
    ASSERT_EQ(1, genotypes.size());

    const auto& sampleGenotype = genotypes.at(0);
    ASSERT_EQ(4, sampleGenotype.data.size());
    EXPECT_EQ("0/1", sampleGenotype.data.at(0).value.get());
    EXPECT_EQ("2", sampleGenotype.data.at(1).value.get());
    EXPECT_EQ("5", sampleGenotype.data.at(2).value.get());
    const auto& acData = sampleGenotype.data.at(3);
    ASSERT_EQ(2, acData.values->size());
    EXPECT_EQ("1", acData.values->at(0));
    EXPECT_EQ("2", acData.values->at(1));

    //    ASSERT_TRUE(sampleGenotype.values.is_initialized());
}
