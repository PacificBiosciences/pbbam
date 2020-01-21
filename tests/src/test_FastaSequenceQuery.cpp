// Author: Derek Barnett

#include <cstddef>

#include <iterator>

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
        FastaSequenceQuery query{fn};
        EXPECT_EQ(1, std::distance(query.begin(), query.end()));
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
        FastaSequenceQuery query{fn};
        EXPECT_EQ(5, std::distance(query.begin(), query.end()));  // 1 from lambda, 4 from chimera
    }
    {
        FastaSequenceQuery query{fn};
        const auto first = query.cbegin();
        const auto& seq = *first;
        EXPECT_EQ("lambda_NEB3011", seq.Name());
    }
}
