// Author: Derek Barnett

#include <pbbam/FastaCache.h>

#include <cctype>

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;

namespace FastaCacheTests {

}  // namespace FastaCacheTests

TEST(FastaCacheTest, can_load_simple)
{
    const std::vector<std::string> expectedNames{
        "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7", "seq8",
    };

    const std::string fn{BAM::PbbamTestsConfig::Data_Dir + "/fastx/simple.fa"};
    const auto cache = BAM::MakeFastaCache(fn);
    EXPECT_EQ(cache->Names(), expectedNames);
    EXPECT_EQ(cache->Subsequence("seq5", 5, 10), "CGTAC");
}

TEST(FastaCacheTest, can_check_sequences)
{
    {
        const std::string fn{BAM::PbbamTestsConfig::Data_Dir + "/fastx/simple.fa"};
        const auto cache = BAM::MakeFastaCache(fn);
        const auto check = cache->Check();
        EXPECT_TRUE(check.first);
        EXPECT_TRUE(check.second.empty());
    }
    {
        const std::string fn{BAM::PbbamTestsConfig::Data_Dir + "/fastx/fasta_cache_check.fa"};
        const auto cache = BAM::MakeFastaCache(fn);
        const auto check = cache->Check();
        EXPECT_FALSE(check.first);
        EXPECT_EQ("bad_seq", check.second);
    }
}

TEST(FastaCacheTest, can_check_sequences_using_callback)
{
    {
        const std::string fn{BAM::PbbamTestsConfig::Data_Dir + "/fastx/simple.fa"};
        const auto cache = BAM::MakeFastaCache(fn);
        const auto check = cache->Check([](const auto& seq) { return seq.Bases().size() == 63; });
        EXPECT_TRUE(check.first);
        EXPECT_TRUE(check.second.empty());
    }
    {
        const std::string fn{BAM::PbbamTestsConfig::Data_Dir + "/fastx/fasta_cache_check.fa"};
        const auto cache = BAM::MakeFastaCache(fn);

        const auto check = cache->Check([](const auto& fasta) {
            const auto& seq = fasta.Bases();
            const auto gcCount = std::count_if(seq.cbegin(), seq.cend(), [](const char base) {
                return (std::toupper(base) == 'C' || std::toupper(base) == 'G');
            });
            return (static_cast<float>(gcCount) / seq.size()) <= 0.5f;
        });

        EXPECT_FALSE(check.first);
        EXPECT_EQ("gc_over_50", check.second);
    }
}
