// Author: Derek Barnett

#include <unistd.h>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/../../src/FileUtils.h>
#include <pbbam/BamFile.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace BamFileTests {

template <typename T>
void CheckFile(const T& input, const size_t expectedCount)
{
    size_t observedCount = 0;
    EntireFileQuery entireFile(input);
    for (const BamRecord& r : entireFile) {
        UNUSED(r);
        ++observedCount;
    }
    EXPECT_EQ(expectedCount, observedCount);
}

}  // namespace BamFileTests

TEST(BamFileTest, NonExistentFileThrows)
{
    EXPECT_THROW(BamFile{"does_not_exist.bam"}, std::runtime_error);
}

TEST(BamFileTest, NonBamFileThrows)
{
    EXPECT_THROW(BamFile{PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa.fai"}, std::runtime_error);
}

TEST(BamFileTest, RelativePathBamOk)
{
    // cache current working directory, then drill down so we can point to
    // BAMs using relative path
    const std::string cwd = internal::FileUtils::CurrentWorkingDirectory();
    ASSERT_EQ(0, chdir(PbbamTestsConfig::Data_Dir.c_str()));
    ASSERT_EQ(0, chdir("relative/a"));

    // BamFile from relative BAM fn
    BamFileTests::CheckFile(BamFile{"../b/test1.bam"}, 3);

    // dataset from relative BAM fn
    BamFileTests::CheckFile(DataSet{"../b/test1.bam"}, 3);

    // dataset from BamFile object (itself from relative BAM fn)
    {
        auto file = BamFile{"../b/test1.bam"};
        BamFileTests::CheckFile(DataSet{file}, 3);
    }

    // restore working directory
    ASSERT_EQ(0, chdir(cwd.c_str()));
}

TEST(BamFileTest, RelativePathXmlOk)
{
    // cache current working directory, then drill down so we can point to
    // BAMs using relative path
    const std::string cwd = internal::FileUtils::CurrentWorkingDirectory();
    ASSERT_EQ(0, chdir(PbbamTestsConfig::Data_Dir.c_str()));

    // dataset from XML containing relative paths
    BamFileTests::CheckFile(DataSet{"relative/relative.xml"}, 9);

    // restore working directory
    ASSERT_EQ(0, chdir(cwd.c_str()));
}

TEST(BamFileTest, RelativePathFofnOk)
{
    // cache current working directory, then drill down so we can point to
    // BAMs using relative path
    const std::string cwd = internal::FileUtils::CurrentWorkingDirectory();
    ASSERT_EQ(0, chdir(PbbamTestsConfig::Data_Dir.c_str()));

    // dataset from FOFN containing relative paths
    BamFileTests::CheckFile(DataSet{"relative/relative.fofn"}, 9);

    // NOTE: doesn't yet support a FOFN containing an XML with relative paths
    //       BamFileTests::CheckFile(DataSet{ "relative/relative2.fofn" }, 60);

    // restore working directory
    ASSERT_EQ(0, chdir(cwd.c_str()));
}

TEST(BamFileTest, TruncatedFileThrowsOk)
{
    EXPECT_THROW(BamFile{PbbamTestsConfig::GeneratedData_Dir + "/truncated.bam"},
                 std::runtime_error);
}
