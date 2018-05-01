// Author: Derek Barnett

#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <boost/any.hpp>

#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/GenomicIntervalQuery.h>
#include <pbbam/Unused.h>
#include <pbbam/ZmwGroupQuery.h>
#include <pbbam/ZmwQuery.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace DataSetQueryTests {

const std::string alignedBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
const std::string aligned2BamFn = PbbamTestsConfig::Data_Dir + "/aligned2.bam";
const std::string alignedCopyBamFn = PbbamTestsConfig::GeneratedData_Dir + "/aligned.bam";
const std::string aligned2CopyBamFn = PbbamTestsConfig::GeneratedData_Dir + "/aligned2.bam";

const std::string group_fofn = PbbamTestsConfig::Generated_Dir + "/group.fofn";
const std::string group_file1 = PbbamTestsConfig::Data_Dir + "/group/test1.bam";
const std::string group_file2 = PbbamTestsConfig::Data_Dir + "/group/test2.bam";
const std::string group_file3 = PbbamTestsConfig::Data_Dir + "/group/test3.bam";

const std::vector<std::string> group_file1_names{
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/24962/0_427"};

const std::vector<std::string> group_file2_names{
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2114_2531",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/4101_5571",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"};

const std::vector<std::string> group_file3_names{
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/45203/0_893",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/45203/0_893",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/46835/3759_4005",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/46835/4052_4686",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/46835/4732_4869",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/47698/9482_9628",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/47698/9675_10333",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/47698/10378_10609",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/49050/48_1132",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/49050/48_1132",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/49194/0_798",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/49194/845_1541",
    "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/49521/0_134"};

static inline bool InGroup(const std::string& name, const std::vector<std::string>& group)
{
    for (const std::string& s : group) {
        if (s == name) return true;
    }
    return false;
}

}  // namespace DataSetQueryTests

TEST(DataSetQueryTest, EntireFileQueryTest)
{
    // single file
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::alignedBamFn);

        DataSet dataset;
        dataset.ExternalResources().Add(bamFile);

        int count = 0;
        EntireFileQuery query(dataset);  // from DataSet object
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(4, count);

        count = 0;
        EntireFileQuery query2(DataSetQueryTests::alignedBamFn);  // from BAM filename
        for (const BamRecord& record : query2) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(4, count);

        count = 0;
        EntireFileQuery query3(bamFile);  // from BamFile object
        for (const BamRecord& record : query3) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(4, count);
    });

    // duplicate file attempt
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::alignedBamFn);

        DataSet dataset;
        dataset.ExternalResources().Add(bamFile);
        dataset.ExternalResources().Add(bamFile);

        int count = 0;
        EntireFileQuery query(dataset);
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(4, count);  // same as single
    });

    // true multi-file dataset
    EXPECT_NO_THROW({
        BamFile file1(DataSetQueryTests::group_file1);  // 1 read
        BamFile file2(DataSetQueryTests::group_file2);  // 4 reads
        BamFile file3(DataSetQueryTests::group_file3);  // 13 reads

        DataSet dataset;
        dataset.ExternalResources().Add(file1);
        dataset.ExternalResources().Add(file2);
        dataset.ExternalResources().Add(file3);

        int count = 0;
        EntireFileQuery query(dataset);
        for (const BamRecord& record : query) {

            // ensure sequential merge of files
            if (count == 0)
                EXPECT_TRUE(DataSetQueryTests::InGroup(record.FullName(),
                                                       DataSetQueryTests::group_file1_names));
            else if (count < 5)
                EXPECT_TRUE(DataSetQueryTests::InGroup(record.FullName(),
                                                       DataSetQueryTests::group_file2_names));
            else
                EXPECT_TRUE(DataSetQueryTests::InGroup(record.FullName(),
                                                       DataSetQueryTests::group_file3_names));

            ++count;
        }
        EXPECT_EQ(18, count);
    });

    // same as above, from FOFN
    EXPECT_NO_THROW({
        int count = 0;

        DataSet dataset(DataSetQueryTests::group_fofn);
        EntireFileQuery query(dataset);
        for (const BamRecord& record : query) {

            // ensure sequential merge of files
            if (count == 0)
                EXPECT_TRUE(DataSetQueryTests::InGroup(record.FullName(),
                                                       DataSetQueryTests::group_file1_names));
            else if (count < 5)
                EXPECT_TRUE(DataSetQueryTests::InGroup(record.FullName(),
                                                       DataSetQueryTests::group_file2_names));
            else
                EXPECT_TRUE(DataSetQueryTests::InGroup(record.FullName(),
                                                       DataSetQueryTests::group_file3_names));

            ++count;
        }
        EXPECT_EQ(18, count);
    });
}

TEST(DataSetQueryTest, GenomicIntervalQueryTest)
{
    const std::string rname = "lambda_NEB3011";

    // single file
    EXPECT_NO_THROW({
        DataSet dataset(DataSetQueryTests::alignedBamFn);  // from BAM filename

        // count records
        int count = 0;
        GenomicInterval interval(rname, 5000, 6000);
        GenomicIntervalQuery query(interval, dataset);
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(2, count);

        // adjust interval and pass back in
        count = 0;
        interval.Start(9000);
        interval.Stop(9500);
        query.Interval(interval);
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(2, count);

        // unknown ref
        count = 0;
        interval.Name("does not exist");
        interval.Start(0);
        interval.Stop(100);
        EXPECT_THROW(query.Interval(interval), std::exception);
        for (const BamRecord& record : query) {  // iteration is still safe, just returns no data
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(0, count);

        // adjust again - make sure we can read a real region after an invalid one
        interval.Name(rname);
        interval.Start(5000);
        interval.Stop(6000);
        query.Interval(interval);
        count = 0;
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(2, count);
    });

    // duplicate file
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::alignedBamFn);

        DataSet dataset;
        dataset.ExternalResources().Add(bamFile);
        dataset.ExternalResources().Add(bamFile);

        // count records & also ensure sorted merge
        int count = 0;
        int prevId = 0;
        int prevPos = 0;

        GenomicInterval interval(rname, 5000, 6000);
        GenomicIntervalQuery query(interval, dataset);
        for (const BamRecord& record : query) {

            EXPECT_TRUE(record.ReferenceId() >= prevId);
            EXPECT_TRUE(record.ReferenceStart() >= prevPos);

            prevId = record.ReferenceId();
            prevPos = record.ReferenceStart();
            ++count;
        }
        EXPECT_EQ(2, count);  // same as single file

        // adjust interval and pass back in
        count = 0;
        interval.Start(9000);
        interval.Stop(10000);
        query.Interval(interval);
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(2, count);  // same as single file

        // unknown ref
        count = 0;
        interval.Name("does not exist");
        interval.Start(0);
        interval.Stop(100);
        EXPECT_THROW(query.Interval(interval), std::exception);
        for (const BamRecord& record : query) {  // iteration is still safe, just returns no data
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(0, count);  // same as single file

        // adjust again - make sure we can read a real region after an invalid one
        interval.Name(rname);
        interval.Start(5000);
        interval.Stop(5300);
        query.Interval(interval);
        count = 0;
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(2, count);  // same as single file
    });

    // multi file BAM (same record content for easy testing, but different filename(ResourceId)
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::alignedBamFn);
        BamFile copyFile(DataSetQueryTests::alignedCopyBamFn);

        DataSet dataset;
        dataset.ExternalResources().Add(bamFile);
        dataset.ExternalResources().Add(copyFile);

        // count records & also ensure sorted merge
        int count = 0;
        int prevId = 0;
        int prevPos = 0;

        GenomicInterval interval(rname, 5000, 6000);
        GenomicIntervalQuery query(interval, dataset);
        for (const BamRecord& record : query) {

            EXPECT_TRUE(record.ReferenceId() >= prevId);
            EXPECT_TRUE(record.ReferenceStart() >= prevPos);

            prevId = record.ReferenceId();
            prevPos = record.ReferenceStart();
            ++count;
        }
        EXPECT_EQ(4, count);  // single file * 2

        // adjust interval and pass back in
        count = 0;
        interval.Start(9000);
        interval.Stop(10000);
        query.Interval(interval);
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(4, count);  // single file * 2

        // unknown ref
        count = 0;
        interval.Name("does not exist");
        interval.Start(0);
        interval.Stop(100);
        EXPECT_THROW(query.Interval(interval), std::exception);
        for (const BamRecord& record : query) {  // iteration is still safe, just returns no data
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(0, count);  // single file * 2

        // adjust again - make sure we can read a real region after an invalid one
        interval.Name(rname);
        interval.Start(5000);
        interval.Stop(5300);
        query.Interval(interval);
        count = 0;
        for (const BamRecord& record : query) {
            UNUSED(record);
            ++count;
        }
        EXPECT_EQ(4, count);  // single file * 2
    });
}

// TODO: implement me
TEST(DataSetQueryTest, QNameQueryTest) { EXPECT_TRUE(true); }

TEST(DataSetQueryTest, ZmwQueryTest)
{
    const std::vector<int32_t> whitelist = {13473, 30983};

    // single file
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::aligned2BamFn);
        ASSERT_TRUE(bamFile.PacBioIndexExists());
        DataSet dataset(bamFile);

        int count = 0;
        ZmwQuery query(whitelist, dataset);
        for (const BamRecord& record : query) {
            const int32_t holeNumber = record.HoleNumber();
            EXPECT_TRUE(holeNumber == 13473 || holeNumber == 30983);
            ++count;
        }
        EXPECT_EQ(4, count);
    });

    // multi-file
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::aligned2BamFn);
        BamFile bamFile2(DataSetQueryTests::aligned2CopyBamFn);
        ASSERT_TRUE(bamFile.PacBioIndexExists());
        ASSERT_TRUE(bamFile2.PacBioIndexExists());

        DataSet dataset;
        dataset.ExternalResources().Add(ExternalResource(bamFile));
        dataset.ExternalResources().Add(ExternalResource(bamFile2));

        int count = 0;
        ZmwQuery query(whitelist, dataset);
        for (const BamRecord& r : query) {
            const auto holeNumber = r.HoleNumber();
            EXPECT_TRUE(holeNumber == 13473 || holeNumber == 30983);
            ++count;
        }
        EXPECT_EQ(8, count);
    });
}

TEST(DataSetQueryTest, ZmwGroupQueryTest)
{
    const std::vector<int32_t> whitelist = {13473, 30983};

    // single-file
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::aligned2BamFn);
        ASSERT_TRUE(bamFile.PacBioIndexExists());
        DataSet dataset(bamFile);

        int count = 0;
        int32_t groupZmw = -1;
        ZmwGroupQuery query(whitelist, dataset);
        for (const std::vector<BamRecord>& group : query) {
            for (const BamRecord& record : group) {
                const auto holeNumber = record.HoleNumber();
                if (groupZmw == -1) groupZmw = holeNumber;
                EXPECT_TRUE(holeNumber == 13473 || holeNumber == 30983);
                EXPECT_EQ(groupZmw, holeNumber);
                ++count;
            }
            groupZmw = -1;
        }
        EXPECT_EQ(4, count);
    });

    // multi-file
    EXPECT_NO_THROW({
        BamFile bamFile(DataSetQueryTests::aligned2BamFn);
        BamFile bamFile2(DataSetQueryTests::aligned2CopyBamFn);
        ASSERT_TRUE(bamFile.PacBioIndexExists());
        ASSERT_TRUE(bamFile2.PacBioIndexExists());

        DataSet dataset;
        dataset.ExternalResources().Add(ExternalResource(bamFile));
        dataset.ExternalResources().Add(ExternalResource(bamFile2));

        int totalCount = 0;
        int numRecordsInGroup = 0;
        int groupCount = 0;
        int32_t groupZmw = -1;
        ZmwGroupQuery query(whitelist, dataset);
        for (const std::vector<BamRecord>& group : query) {
            for (const BamRecord& record : group) {
                const auto holeNumber = record.HoleNumber();
                ++numRecordsInGroup;
                if (groupZmw == -1) groupZmw = holeNumber;
                EXPECT_TRUE(holeNumber == 13473 || holeNumber == 30983);
                EXPECT_EQ(groupZmw, holeNumber);
                ++totalCount;
            }
            if (groupCount == 0)
                EXPECT_EQ(4, numRecordsInGroup);
            else if (groupCount == 1)
                EXPECT_EQ(4, numRecordsInGroup);
            else
                EXPECT_TRUE(false);  // should not get here
            numRecordsInGroup = 0;
            ++groupCount;
            groupZmw = -1;
        }
        EXPECT_EQ(8, totalCount);
    });
}
