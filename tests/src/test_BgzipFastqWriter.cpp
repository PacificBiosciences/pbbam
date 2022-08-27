#include <pbbam/BgzipFastqWriter.h>

#include <string>

#include <gtest/gtest.h>

#include <pbbam/FastqReader.h>
#include <pbbam/FastqSequence.h>
#include <pbbam/FormatUtils.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_BgzipFastqWriter, writes_bgzf_fastq_data_to_file)
{
    const std::string fn{PbbamTestsConfig::GeneratedData_Dir + "/bgzf_fastq_out.fq.gz"};

    const std::vector<FastqSequence> sequences{
        FastqSequence{"seq1", "ACGT", Data::QualityValues{"zzzz"}},
        FastqSequence{"seq2", "GATTACA", Data::QualityValues{"~~~~~~~"}},
        FastqSequence{"seq3", "CCCC", Data::QualityValues{"$$$$"}}};

    {
        BgzipFastqWriter writer{fn};
        for (const auto& seq : sequences) {
            writer.Write(seq);
        }
    }
    EXPECT_EQ(HtslibCompression::BGZIP, FormatUtils::CompressionType(fn));

    const auto observed = FastqReader::ReadAll(fn);
    EXPECT_TRUE(std::equal(sequences.cbegin(), sequences.cend(), observed.cbegin()));
}
