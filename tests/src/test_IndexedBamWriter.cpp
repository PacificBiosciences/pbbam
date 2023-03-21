#include <pbbam/IndexedBamWriter.h>

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiBuilder.h>
#include <pbbam/PbiRawData.h>

#include "PbbamTestData.h"

TEST(BAM_IndexedBamWriter, writes_valid_bam_and_pbi_on_success)
{
    using namespace PacBio::BAM;

    const std::string inBam = PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam";
    const std::string outBam = PbbamTestsConfig::GeneratedData_Dir + "/ibw.bam";
    const std::string outPbi = PbbamTestsConfig::GeneratedData_Dir + "/ibw.bam.pbi";

    const BamFile file{inBam};
    const auto& header = file.Header();

    const std::vector<std::string> expectedQNames{
        "ArminsFakeMovie/100000/2659_3025", "ArminsFakeMovie/100000/3116_3628",
        "ArminsFakeMovie/100000/3722_4267", "ArminsFakeMovie/100000/4356_4864",
        "ArminsFakeMovie/100000/4960_5477", "ArminsFakeMovie/100000/5571_6087",
        "ArminsFakeMovie/100000/6199_6719", "ArminsFakeMovie/100000/6812_7034",
        "ArminsFakeMovie/200000/2659_3025", "ArminsFakeMovie/200000/3116_3628",
        "ArminsFakeMovie/200000/3722_4267", "ArminsFakeMovie/200000/4356_4864",
        "ArminsFakeMovie/200000/4960_5477", "ArminsFakeMovie/200000/5571_6087",
        "ArminsFakeMovie/200000/6199_6719", "ArminsFakeMovie/200000/6812_7034",
        "ArminsFakeMovie/300000/2659_3025", "ArminsFakeMovie/300000/3116_3628",
        "ArminsFakeMovie/300000/3722_4267", "ArminsFakeMovie/300000/4356_4864",
        "ArminsFakeMovie/300000/4960_5477", "ArminsFakeMovie/300000/5571_6087",
        "ArminsFakeMovie/300000/6199_6719", "ArminsFakeMovie/300000/6812_7034"};

    {  // copy file & generate index

        BamReader reader{file};
        IndexedBamWriter writer{outBam, header};
        BamRecord b;
        while (reader.GetNext(b)) {
            writer.Write(b);
        }
    }

    {  // sequential read of new BAM

        BamReader reader{outBam};
        BamRecord b;
        for (size_t i = 0; i < 24; ++i) {
            reader.GetNext(b);
            EXPECT_EQ(expectedQNames.at(i), b.FullName());
        }
    }

    {  // check random access in new BAM, using companion PBI

        const PbiRawData idx{outPbi};
        const auto& offsets = idx.BasicData().fileOffset_;

        BamReader reader{outBam};
        BamRecord b;
        for (int i = 23; i >= 0; --i) {
            reader.VirtualSeek(offsets.at(i));
            reader.GetNext(b);
            EXPECT_EQ(expectedQNames.at(i), b.FullName());
        }
    }
}

TEST(BAM_IndexedBamWriter, can_handle_long_reads_spanning_bgzf_blocks)
{
    using namespace PacBio::BAM;

    const std::string inBamFn = PbbamTestsConfig::Data_Dir + "/long_reads.bam";
    const std::string outBamFn = PbbamTestsConfig::GeneratedData_Dir + "/long_reads.copy.bam";
    const std::string outPbiFn = PbbamTestsConfig::GeneratedData_Dir + "/long_reads.copy.bam.pbi";

    // copy file, writing inline PBI index
    {
        BamFile file{inBamFn};
        IndexedBamWriter writer{outBamFn, file.Header()};
        EntireFileQuery query{file};
        for (const auto& b : query) {
            writer.Write(b);
        }
    }

    {
        const PbiRawData idx{outPbiFn};
        const auto& offsets = idx.BasicData().fileOffset_;

        BamReader reader{outBamFn};
        BamRecord b;
        for (int i = 0; i < 100; ++i) {
            reader.VirtualSeek(offsets.at(i));
            reader.GetNext(b);
        }
    }

    remove(outBamFn.c_str());
    remove(outPbiFn.c_str());
}

TEST(BAM_IndexedBamWriter, removes_gzi_file_for_bam_with_no_records)
{
    const std::string inBamFn = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/long_reads.bam";
    const std::string outBamFn =
        PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/long_reads.copy.bam";

    PacBio::BAM::BamFile file{inBamFn};

    // GZI file removed for normal, non-empty BAM (header + records)
    {
        PacBio::BAM::IndexedBamWriter writer{outBamFn, file.Header()};
        PacBio::BAM::EntireFileQuery query{file};
        for (const auto& record : query) {
            writer.Write(record);
        }
    }
    EXPECT_TRUE(std::filesystem::exists(outBamFn));
    EXPECT_TRUE(std::filesystem::exists(outBamFn + ".pbi"));
    EXPECT_FALSE(std::filesystem::exists(outBamFn + ".gzi"));
    std::filesystem::remove(outBamFn);
    std::filesystem::remove(outBamFn + ".pbi");

    // GZI file removed for empty BAM (header-only)
    {
        PacBio::BAM::IndexedBamWriter writer{outBamFn, file.Header()};
    }
    EXPECT_TRUE(std::filesystem::exists(outBamFn));
    EXPECT_TRUE(std::filesystem::exists(outBamFn + ".pbi"));
    EXPECT_FALSE(std::filesystem::exists(outBamFn + ".gzi"));
    std::filesystem::remove(outBamFn);
    std::filesystem::remove(outBamFn + ".pbi");
}
