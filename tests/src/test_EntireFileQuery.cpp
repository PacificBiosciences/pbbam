// Author: Derek Barnett

#include <pbbam/EntireFileQuery.h>

#include <iterator>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/BamWriter.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace EntireFileQueryTests {

const std::string inputBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";

}  // namespace EntireFileQueryTests

TEST(EntireFileQueryTest, CountRecords)
{
    const BamFile bamFile{EntireFileQueryTests::inputBamFn};
    EntireFileQuery entireFile{bamFile};
    EXPECT_EQ(4, std::distance(entireFile.begin(), entireFile.end()));
}

TEST(BamRecordTest, HandlesDeletionOK)
{
    // this file raised no error in Debug mode, but segfaulted when
    // trying to access the aligned qualities in Release mode
    const BamFile bamFile{PbbamTestsConfig::Data_Dir + "/segfault.bam"};

    int count = 0;
    EntireFileQuery entireFile{bamFile};
    for (const BamRecord& record : entireFile) {

        const auto rawQualities = record.Qualities(Orientation::GENOMIC, false);
        const auto alignedQualities = record.Qualities(Orientation::GENOMIC, true);

        const std::string rawExpected{
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
            "IIIIIIIIIIIII"};

        // 1=1D98=
        const std::string alignedExpected{
            "I!"
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
            "IIIIIIIIIIII"};

        EXPECT_EQ(rawExpected, rawQualities.Fastq());
        EXPECT_EQ(alignedExpected, alignedQualities.Fastq());

        ++count;
    }

    EXPECT_EQ(1, count);
}

TEST(BamRecordTest, ReferenceName)
{
    {  // check reference name of first record
        const BamFile bamFile{PbbamTestsConfig::Data_Dir + "/aligned.bam"};
        EntireFileQuery records{bamFile};
        auto firstIter = records.begin();
        auto& firstRecord = *firstIter;
        ASSERT_TRUE(firstRecord.IsMapped());
        EXPECT_EQ("lambda_NEB3011", firstRecord.ReferenceName());
    }

    {  // unmapped records have no reference name, should throw
        const BamFile bamFile{PbbamTestsConfig::Data_Dir + "/unmap1.bam"};
        EntireFileQuery records{bamFile};
        auto firstIter = records.begin();
        auto& firstRecord = *firstIter;
        ASSERT_FALSE(firstRecord.IsMapped());
        EXPECT_THROW(firstRecord.ReferenceName(), std::runtime_error);
    }
}
