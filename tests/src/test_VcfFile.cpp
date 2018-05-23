// Author: Derek Barnett

#include <gtest/gtest.h>

#include <pbbam/vcf/VcfFile.h>
#include <pbbam/vcf/VcfFormat.h>

#include "PbbamTestData.h"

using VcfFile = PacBio::VCF::VcfFile;
using VcfFormat = PacBio::VCF::VcfFormat;

namespace VcfFileTests {

static const std::string BasicHeaderText{
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

static const std::string VcfFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                               "/vcf/structural_variants.vcf"};

}  // namespace VcfFileTests

TEST(VCF_File, initializes_header_from_input_file)
{
    const VcfFile file{VcfFileTests::VcfFn};
    const auto hdrText = VcfFormat::FormattedHeader(file.Header());

    EXPECT_EQ(VcfFileTests::BasicHeaderText, hdrText);
}
