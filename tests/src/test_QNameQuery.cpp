// Author: Yuan Li

#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/QNameQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace QNameQueryTests {

static const string dataDir = PbbamTestsConfig::Data_Dir + "/group/";
static const string test1fn = string(dataDir) + "test1.bam";
static const string test2fn = string(dataDir) + "test2.bam";
static const string test3fn = string(dataDir) + "test3.bam";

static void TestQNameQuery(const string& fn, const vector<int>& expected)
{
    EXPECT_NO_THROW({
        vector<int> counts;
        QNameQuery qQuery(fn);
        for (const vector<BamRecord>& records : qQuery)
            counts.push_back(records.size());
        EXPECT_EQ(expected, counts);
    });
}

static void TestNoneConstQNameQuery(const string& fn, const vector<int>& expected)
{
    EXPECT_NO_THROW({
        vector<int> counts;
        QNameQuery qQuery(fn);
        for (vector<BamRecord>& records : qQuery)
            counts.push_back(records.size());
        EXPECT_EQ(expected, counts);
    });
}

}  // namespace QNameQueryTests

TEST(QNameQueryTest, CountQSizes)
{
    // test case 1 has exactly one bamRecord.
    string fn = QNameQueryTests::test1fn;
    vector<int> expected({1});
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
