// Author: Derek Barnett

#include <pbbam/FaiIndex.h>

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace FaiIndexTests {

const std::string simpleFastaFn{PbbamTestsConfig::Data_Dir + "/fastx/simple.fa"};
const std::string simpleFastaFaiFn{PbbamTestsConfig::Data_Dir + "/fastx/simple.fa.fai"};
const std::string simpleFastqFn{PbbamTestsConfig::Data_Dir + "/fastx/simple.fq"};
const std::string simpleFastqFaiFn{PbbamTestsConfig::Data_Dir + "/fastx/simple.fq.fai"};

}  // namespace FaiIndexTests

TEST(FaiIndexTest, LoadsFromFastaFaiFile)
{
    const std::vector<std::string> expectedNames{"seq1", "seq2", "seq3", "seq4",
                                                 "seq5", "seq6", "seq7", "seq8"};
    const FaiEntry expectedEntry{63, 146, 63, 64};

    const FaiIndex index{FaiIndexTests::simpleFastaFaiFn};
    const auto& names = index.Names();
    ASSERT_EQ(8, names.size());
    EXPECT_TRUE(std::equal(expectedNames.cbegin(), expectedNames.cend(), names.cbegin()));
    EXPECT_EQ(expectedEntry, index.Entry("seq3"));
}

TEST(FaiIndexTest, LoadsFromFastqFaiFile)
{
    const std::vector<std::string> expectedNames{"seq1", "seq2", "seq3", "seq4",
                                                 "seq5", "seq6", "seq7", "seq8"};
    const FaiEntry expectedEntry{63, 278, 63, 64, 344};

    const FaiIndex index{FaiIndexTests::simpleFastqFaiFn};
    const auto& names = index.Names();
    ASSERT_EQ(8, names.size());
    EXPECT_TRUE(std::equal(expectedNames.cbegin(), expectedNames.cend(), names.cbegin()));
    EXPECT_EQ(expectedEntry, index.Entry("seq3"));
}

TEST(FaiIndexTest, SaveFastaIndexToStream)
{
    // clang-format off
    const std::string expected
    {
        "seq1\t63\t6\t63\t64\n"
        "seq2\t63\t76\t63\t64\n"
        "seq3\t63\t146\t63\t64\n"
        "seq4\t63\t216\t63\t64\n"
        "seq5\t63\t286\t63\t64\n"
        "seq6\t63\t356\t63\t64\n"
        "seq7\t63\t426\t63\t64\n"
        "seq8\t63\t496\t63\t64\n"
    };
    // clang-format on

    FaiIndex index;
    index.Add("seq1", {63, 6, 63, 64});
    index.Add("seq2", {63, 76, 63, 64});
    index.Add("seq3", {63, 146, 63, 64});
    index.Add("seq4", {63, 216, 63, 64});
    index.Add("seq5", {63, 286, 63, 64});
    index.Add("seq6", {63, 356, 63, 64});
    index.Add("seq7", {63, 426, 63, 64});
    index.Add("seq8", {63, 496, 63, 64});

    std::ostringstream out;
    index.Save(out);
    EXPECT_EQ(expected, out.str());
}

TEST(FaiIndexTest, SaveFastqIndexToStream)
{
    // clang-format off
    const std::string expected
    {
        "seq1\t63\t6\t63\t64\t72\n"
        "seq2\t63\t142\t63\t64\t208\n"
        "seq3\t63\t278\t63\t64\t344\n"
        "seq4\t63\t414\t63\t64\t480\n"
        "seq5\t63\t550\t63\t64\t616\n"
        "seq6\t63\t686\t63\t64\t752\n"
        "seq7\t63\t822\t63\t64\t888\n"
        "seq8\t63\t958\t63\t64\t1024\n"
    };
    // clang-format on

    FaiIndex index;
    index.Add("seq1", {63, 6, 63, 64, 72});
    index.Add("seq2", {63, 142, 63, 64, 208});
    index.Add("seq3", {63, 278, 63, 64, 344});
    index.Add("seq4", {63, 414, 63, 64, 480});
    index.Add("seq5", {63, 550, 63, 64, 616});
    index.Add("seq6", {63, 686, 63, 64, 752});
    index.Add("seq7", {63, 822, 63, 64, 888});
    index.Add("seq8", {63, 958, 63, 64, 1024});

    std::ostringstream out;
    index.Save(out);
    EXPECT_EQ(expected, out.str());
}

TEST(FaiIndexTest, throws_on_missing_fai_file)
{
    EXPECT_THROW(FaiIndex{"does_not_exist.fai"}, std::runtime_error);
}
