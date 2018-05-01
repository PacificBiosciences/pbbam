// Author: Yuan Li

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/QNameQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace QNameQueryTests {

static const std::string dataDir = PbbamTestsConfig::Data_Dir + "/group/";
static const std::string test1fn = std::string(dataDir) + "test1.bam";
static const std::string test2fn = std::string(dataDir) + "test2.bam";
static const std::string test3fn = std::string(dataDir) + "test3.bam";

static void TestQNameQuery(const std::string& fn, const std::vector<int>& expected)
{
    EXPECT_NO_THROW({
        std::vector<int> counts;
        QNameQuery qQuery(fn);
        for (const std::vector<BamRecord>& records : qQuery)
            counts.push_back(records.size());
        EXPECT_EQ(expected, counts);
    });
}

static void TestNoneConstQNameQuery(const std::string& fn, const std::vector<int>& expected)
{
    EXPECT_NO_THROW({
        std::vector<int> counts;
        QNameQuery qQuery(fn);
        for (std::vector<BamRecord>& records : qQuery)
            counts.push_back(records.size());
        EXPECT_EQ(expected, counts);
    });
}

}  // namespace QNameQueryTests

TEST(QNameQueryTest, CountQSizes)
{
    // test case 1 has exactly one bamRecord.
    std::string fn = QNameQueryTests::test1fn;
    std::vector<int> expected({1});
    QNameQueryTests::TestQNameQuery(fn, expected);
    QNameQueryTests::TestNoneConstQNameQuery(fn, expected);

    // test case 2 has bamRecords of four subreads.
    fn = QNameQueryTests::test2fn;
    expected = {1, 1, 1, 1};
    QNameQueryTests::TestQNameQuery(fn, expected);
    QNameQueryTests::TestNoneConstQNameQuery(fn, expected);

    fn = QNameQueryTests::test3fn;
    expected = {2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1};
    QNameQueryTests::TestQNameQuery(fn, expected);
    QNameQueryTests::TestNoneConstQNameQuery(fn, expected);
}
