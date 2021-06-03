#include <pbbam/BamFile.h>

#include <cstddef>
#include <cstdlib>

#include <iterator>
#include <stdexcept>

#include <unistd.h>

#include <gtest/gtest.h>

#include <pbbam/../../src/FileUtils.h>
#include <pbbam/EntireFileQuery.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace BamFileTests {

template <typename T>
void CheckFile(const T& input, const size_t expectedCount)
{
    EntireFileQuery entireFile{input};
    EXPECT_EQ(expectedCount, std::distance(entireFile.begin(), entireFile.end()));
}

}  // namespace BamFileTests

TEST(BAM_BamFile, throws_on_non_existent_file)
{
    EXPECT_THROW(BamFile{"does_not_exist.bam"}, std::runtime_error);
}

TEST(BAM_BamFile, throws_on_wrong_file_format)
{
    EXPECT_THROW(BamFile{PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa.fai"}, std::runtime_error);
}

TEST(BAM_BamFile, throws_on_truncated_file)
{
    EXPECT_THROW(BamFile{PbbamTestsConfig::GeneratedData_Dir + "/truncated.bam"},
                 std::runtime_error);
}

TEST(BAM_BamFile, can_load_from_relative_path_bam)
{
    // cache current working directory, then drill down so we can point to
    // BAMs using relative path
    const std::string cwd{FileUtils::CurrentWorkingDirectory()};
    ASSERT_EQ(0, chdir(PbbamTestsConfig::Data_Dir.c_str()));
    ASSERT_EQ(0, chdir("relative/a"));

    // BamFile from relative BAM fn
    BamFileTests::CheckFile(BamFile{"../b/test1.bam"}, 3);

    // dataset from relative BAM fn
    BamFileTests::CheckFile(DataSet{"../b/test1.bam"}, 3);

    // dataset from BamFile object (itself from relative BAM fn)
    {
        const BamFile file{"../b/test1.bam"};
        BamFileTests::CheckFile(DataSet{file}, 3);
    }

    // restore working directory
    ASSERT_EQ(0, chdir(cwd.c_str()));
}

TEST(BAM_BamFile, can_load_from_relative_path_dataset)
{
    // cache current working directory, then drill down so we can point to
    // BAMs using relative path
    const std::string cwd{FileUtils::CurrentWorkingDirectory()};
    ASSERT_EQ(0, chdir(PbbamTestsConfig::Data_Dir.c_str()));

    // dataset from XML containing relative paths
    BamFileTests::CheckFile(DataSet{"relative/relative.xml"}, 9);

    // restore working directory
    ASSERT_EQ(0, chdir(cwd.c_str()));
}

TEST(BAM_BamFile, can_load_from_relative_path_fofn)
{
    // cache current working directory, then drill down so we can point to
    // BAMs using relative path
    const std::string cwd{FileUtils::CurrentWorkingDirectory()};
    ASSERT_EQ(0, chdir(PbbamTestsConfig::Data_Dir.c_str()));

    // dataset from FOFN containing relative paths
    BamFileTests::CheckFile(DataSet{"relative/relative.fofn"}, 9);

    // NOTE: doesn't yet support a FOFN containing an XML with relative paths
    //       BamFileTests::CheckFile(DataSet{ "relative/relative2.fofn" }, 60);

    // restore working directory
    ASSERT_EQ(0, chdir(cwd.c_str()));
}
