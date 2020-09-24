// Author: Derek Barnett

#include <pbbam/FastaSequenceQuery.h>

#include <iterator>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_FastaSequenceQuery, can_read_from_fasta_file)
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

TEST(BAM_FastaSequenceQuery, can_read_from_dataset)
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
