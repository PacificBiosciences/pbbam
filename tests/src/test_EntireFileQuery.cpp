// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace EntireFileQueryTests {

const std::string inputBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";

}  // namespace EntireFileQueryTests

TEST(EntireFileQueryTest, CountRecords)
{
    EXPECT_NO_THROW({
        BamFile bamFile(EntireFileQueryTests::inputBamFn);
        int count = 0;
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile) {
            UNUSED(record);
            ++count;
        }

        EXPECT_EQ(4, count);
    });
}

TEST(EntireFileQueryTest, NonConstBamRecord)
{
    EXPECT_NO_THROW({
        BamFile bamFile(EntireFileQueryTests::inputBamFn);
        int count = 0;
        EntireFileQuery entireFile(bamFile);
        for (BamRecord& record : entireFile) {
            UNUSED(record);
            ++count;
        }

        EXPECT_EQ(4, count);
    });
}

TEST(BamRecordTest, HandlesDeletionOK)
{
    // this file raised no error in Debug mode, but segfaulted when
    // trying to access the aligned qualities in Release mode

    const std::string problemBamFn = PbbamTestsConfig::Data_Dir + "/segfault.bam";
    BamFile bamFile(problemBamFn);
    int count = 0;
    EntireFileQuery entireFile(bamFile);
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
        const std::string exampleBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";
        BamFile bamFile(exampleBam);
        EntireFileQuery records(bamFile);
        auto firstIter = records.begin();
        auto& firstRecord = *firstIter;
        ASSERT_TRUE(firstRecord.IsMapped());
        EXPECT_EQ("lambda_NEB3011", firstRecord.ReferenceName());
    }

    {  // unmapped records have no reference name, should throw
        const std::string exampleBam = PbbamTestsConfig::Data_Dir + "/unmap1.bam";
        BamFile bamFile(exampleBam);
        EntireFileQuery records(bamFile);
        auto firstIter = records.begin();
        auto& firstRecord = *firstIter;
        ASSERT_FALSE(firstRecord.IsMapped());
        EXPECT_THROW(firstRecord.ReferenceName(), std::runtime_error);
    }
}
