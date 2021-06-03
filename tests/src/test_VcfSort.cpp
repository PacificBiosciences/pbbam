#include <pbbam/vcf/VcfSort.h>

#include <cstdio>

#include <gtest/gtest.h>

#include <pbbam/vcf/VcfQuery.h>

#include "PbbamTestData.h"

using VcfFile = PacBio::VCF::VcfFile;
using VcfQuery = PacBio::VCF::VcfQuery;

namespace VcfSortTests {

const std::string inputFn = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/vcf/unsorted.vcf";
const std::string outputFn = PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/sorted.vcf";

}  // namespace VcfSortTests

TEST(VCF_VcfSort, sorts_input_file)
{
    const VcfFile file{VcfSortTests::inputFn};
    PacBio::VCF::SortFile(file, VcfSortTests::outputFn);

    const std::vector<std::string> expectedIds{"variant0", "variant5", "variant1",
                                               "variant3", "variant4", "variant2"};

    size_t i = 0;
    VcfQuery query{VcfSortTests::outputFn};
    for (const auto& var : query) {
        EXPECT_EQ(expectedIds.at(i), var.Id());
        ++i;
    }

    // remove temp file
    remove(VcfSortTests::outputFn.c_str());
}
