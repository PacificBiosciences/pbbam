// Author: Derek Barnett

#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <pbbam/vcf/VcfFile.h>
#include <pbbam/vcf/VcfReader.h>

#include "PbbamTestData.h"

using VcfFile = PacBio::VCF::VcfFile;
using VcfReader = PacBio::VCF::VcfReader;
using VcfVariant = PacBio::VCF::VcfVariant;

namespace VcfReaderTests {

static const std::vector<std::string> ExpectedIds{
    "pbsv.INS.1",  "pbsv.DEL.2",  "pbsv.INS.3",  "pbsv.INS.4",  "pbsv.DEL.5",  "pbsv.DEL.6",
    "pbsv.DEL.7",  "pbsv.INS.8",  "pbsv.INS.9",  "pbsv.INS.10", "pbsv.INS.11", "pbsv.INS.12",
    "pbsv.INS.13", "pbsv.INS.14", "pbsv.INS.15", "pbsv.INS.16", "pbsv.INS.17", "pbsv.INS.18",
    "pbsv.INS.19", "pbsv.DEL.20", "pbsv.INS.21"};

static const std::string VcfFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                               "/vcf/structural_variants.vcf"};

}  // namespace VcfReaderTests

TEST(VCF_Reader, can_fetch_variants_from_vcf_filename)
{
    size_t i = 0;
    VcfReader rdr{VcfReaderTests::VcfFn};
    VcfVariant var;
    while (rdr.GetNext(var)) {
        EXPECT_EQ(VcfReaderTests::ExpectedIds.at(i), var.Id());
        ++i;
    }
}

TEST(VCF_Reader, can_fetch_variants_from_vcf_file_object)
{
    const VcfFile file{VcfReaderTests::VcfFn};

    size_t i = 0;
    VcfReader rdr{file};
    VcfVariant var;
    while (rdr.GetNext(var)) {
        EXPECT_EQ(VcfReaderTests::ExpectedIds.at(i), var.Id());
        ++i;
    }
}
