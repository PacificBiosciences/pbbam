// Author: Derek Barnett

#include <cstddef>

#include <tuple>

#include <gtest/gtest.h>

#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/FastaWriter.h>
#include <boost/algorithm/string.hpp>

#include "FastxTests.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace FastaSequenceQueryTests {
}  // namespace FastaSequenceQueryTests

TEST(FastaSequenceQueryTest, can_read_from_fasta_file)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";

    {
        size_t count = 0;
        FastaSequenceQuery query{fn};
        for (const auto& seq : query) {
            std::ignore = seq;
            ++count;
        }
        EXPECT_EQ(1, count);
    }

    {
        FastaSequenceQuery query{fn};
        const auto first = query.cbegin();
        const auto& seq = *first;
        EXPECT_EQ("lambda_NEB3011", seq.Name());
    }
}

TEST(FastaSequenceQueryTest, can_read_from_dataset)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/referenceset.xml";

    {
        size_t count = 0;
        FastaSequenceQuery query{fn};
        for (const auto& seq : query) {
            std::ignore = seq;
            ++count;
        }
        EXPECT_EQ(5, count);  // 1 from lambda, 4 from chimera
    }
    {
        FastaSequenceQuery query{fn};
        const auto first = query.cbegin();
        const auto& seq = *first;
        EXPECT_EQ("lambda_NEB3011", seq.Name());
    }
}
