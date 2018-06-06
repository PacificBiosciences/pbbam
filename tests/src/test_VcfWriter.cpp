// Author: Derek Barnett

#include <cstdio>

#include <gtest/gtest.h>
#include <pbbam/vcf/VcfFile.h>
#include <pbbam/vcf/VcfFormat.h>
#include <pbbam/vcf/VcfQuery.h>
#include <pbbam/vcf/VcfWriter.h>

#include "PbbamTestData.h"

using VcfFile = PacBio::VCF::VcfFile;
using VcfFormat = PacBio::VCF::VcfFormat;
using VcfQuery = PacBio::VCF::VcfQuery;
using VcfWriter = PacBio::VCF::VcfWriter;

namespace VcfWriterTests {

static const std::string VcfFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                               "/vcf/structural_variants.vcf"};

}  // namespace VcfWriterTests

TEST(VCF_Writer, correctly_copies_vcf_file)
{
    const std::string intitialFn{VcfWriterTests::VcfFn};
    const std::string newFn{PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/temp.vcf"};

    const VcfFile initialFile{VcfWriterTests::VcfFn};

    const std::string expectedHeaderText = VcfFormat::FormattedHeader(initialFile.Header());
    std::vector<std::string> expectedVariantsText;

    {  // store contents of intitial file & write to a new file
        VcfWriter writer{newFn, initialFile.Header()};
        VcfQuery query{initialFile};
        for (const auto& var : query) {
            expectedVariantsText.push_back(VcfFormat::FormattedVariant(var));
            writer.Write(var);
        }
    }
    {  // read new file & compare against original

        const VcfFile newFile{newFn};
        EXPECT_EQ(expectedHeaderText, VcfFormat::FormattedHeader(newFile.Header()));

        size_t i = 0;
        for (const auto& var : VcfQuery{newFile}) {
            EXPECT_EQ(expectedVariantsText.at(i), VcfFormat::FormattedVariant(var));
            ++i;
        }
    }
    ::remove(newFn.c_str());
}
