#include <pbbam/BgzipFastaWriter.h>

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FormatUtils.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_BgzipFastaWriter, writes_bgzf_fasta_data_to_file)
{
    const std::string fn{PbbamTestsConfig::GeneratedData_Dir + "/bgzf_fasta_out.fa.gz"};

    const std::vector<FastaSequence> sequences{FastaSequence{"seq1", "ACGT"},
                                               FastaSequence{"seq2", "GATTACA"},
                                               FastaSequence{"seq3", "CCCC"}};

    {
        BgzipFastaWriter writer{fn};
        for (const auto& seq : sequences) {
            writer.Write(seq);
        }
    }
    EXPECT_EQ(HtslibCompression::BGZIP, FormatUtils::CompressionType(fn));

    const auto observed = FastaReader::ReadAll(fn);
    EXPECT_TRUE(std::equal(sequences.cbegin(), sequences.cend(), observed.cbegin()));
}
