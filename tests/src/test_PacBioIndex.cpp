// Author: Derek Barnett

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiBuilder.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace PacBioIndexTests {

const std::string test2BamFn = PbbamTestsConfig::Data_Dir + "/aligned2.bam";
const std::string phi29BamFn = PbbamTestsConfig::Data_Dir + "/phi29.bam";

static PbiRawData Test2Bam_CoreIndexData()
{
    PbiRawData rawData;
    rawData.Version(PbiFile::Version_3_0_1);
    rawData.FileSections(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::REFERENCE);
    rawData.NumReads(10);

    PbiRawBasicData& basicData = rawData.BasicData();
    basicData.rgId_ = {-1197849594, -1197849594, -1197849594, -1197849594, -1197849594,
                       -1197849594, -1197849594, -1197849594, -1197849594, -1197849594};
    basicData.qStart_ = {48, 387, 0, 9936, 10232, 7468, 5557, 7285, 426, 7064};
    basicData.qEnd_ = {1132, 1134, 344, 10187, 10394, 8906, 7235, 8657, 1045, 7421};
    basicData.holeNumber_ = {49050, 32328, 32328, 6469, 6469, 30983, 13473, 13473, 19915, 30983};
    basicData.readQual_ = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6};
    basicData.ctxtFlag_ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    basicData.fileOffset_ = {33816576, 33825163, 33831333, 33834264, 33836542,
                             33838065, 33849818, 33863499, 33874621, 1392836608};

    PbiRawMappedData& mappedData = rawData.MappedData();
    mappedData.tId_ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    mappedData.tStart_ = {0, 302, 675, 2170, 2203, 3572, 4506, 4507, 4592, 4669};
    mappedData.tEnd_ = {471, 1019, 1026, 2397, 2326, 5015, 6125, 5850, 5203, 5011};
    mappedData.aStart_ = {653, 395, 1, 9960, 10271, 7468, 5574, 7285, 441, 7075};
    mappedData.aEnd_ = {1129, 1134, 344, 10185, 10394, 8906, 7235, 8647, 1040, 7418};
    mappedData.revStrand_ = {0, 1, 0, 1, 0, 1, 1, 0, 1, 0};
    mappedData.nM_ = {460, 704, 339, 216, 118, 1394, 1581, 1313, 583, 333};
    mappedData.nMM_ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    mappedData.mapQV_ = {254, 254, 254, 254, 254, 254, 254, 254, 254, 254};

    PbiRawReferenceData& referenceData = rawData.ReferenceData();
    referenceData.entries_ = {PbiReferenceEntry{0, 0, 10},
                              PbiReferenceEntry{4294967295, 4294967295, 4294967295}};

    return rawData;
}

// NOTE: We have 2 different sets of offsets because the copied, new file differs in size than the existing one.
//
//       Unsure which combination of write parameters were used on the original. Things like thread count,
//       compression level, etc. can effect compression ratio, BGZF block sizes, etc. even though the BAM record
//       content itself is equal. So we'll just track these index values separately, for now at least.
//
static PbiRawData Test2Bam_ExistingIndex()
{
    PbiRawData index = Test2Bam_CoreIndexData();
    index.BasicData().fileOffset_ = {33816576, 33825163, 33831333, 33834264, 33836542,
                                     33838065, 33849818, 33863499, 33874621, 1392836608};
    return index;
}

static PbiRawData Test2Bam_NewIndex()
{
    PbiRawData index = Test2Bam_CoreIndexData();
    index.BasicData().fileOffset_ = {33816576,  236126208, 391315456,  469106688,  537067520,
                                     587792384, 867303424, 1182793728, 1449787392, 1582628864};
    return index;
}

static void ExpectRawIndicesEqual(const PbiRawData& expected, const PbiRawData& actual)
{
    // header data
    EXPECT_EQ(expected.Version(), actual.Version());
    EXPECT_EQ(expected.FileSections(), actual.FileSections());
    EXPECT_EQ(expected.NumReads(), actual.NumReads());

    // subread data
    const PbiRawBasicData& e = expected.BasicData();
    const PbiRawBasicData& a = actual.BasicData();
    EXPECT_EQ(e.rgId_, a.rgId_);
    EXPECT_EQ(e.qStart_, a.qStart_);
    EXPECT_EQ(e.qEnd_, a.qEnd_);
    EXPECT_EQ(e.holeNumber_, a.holeNumber_);
    EXPECT_EQ(e.readQual_, a.readQual_);
    EXPECT_EQ(e.ctxtFlag_, a.ctxtFlag_);
    EXPECT_EQ(e.fileOffset_, a.fileOffset_);

    // mapped data
    EXPECT_EQ(expected.HasMappedData(), actual.HasMappedData());
    if (expected.HasMappedData() && actual.HasMappedData()) {
        const PbiRawMappedData& e2 = expected.MappedData();
        const PbiRawMappedData& a2 = actual.MappedData();
        EXPECT_EQ(e2.tId_, a2.tId_);
        EXPECT_EQ(e2.tStart_, a2.tStart_);
        EXPECT_EQ(e2.tEnd_, a2.tEnd_);
        EXPECT_EQ(e2.aStart_, a2.aStart_);
        EXPECT_EQ(e2.aEnd_, a2.aEnd_);
        EXPECT_EQ(e2.revStrand_, a2.revStrand_);
        EXPECT_EQ(e2.nM_, a2.nM_);
        EXPECT_EQ(e2.nMM_, a2.nMM_);
        EXPECT_EQ(e2.mapQV_, a2.mapQV_);
    }

    // reference data
    EXPECT_EQ(expected.HasReferenceData(), actual.HasReferenceData());
    if (expected.HasReferenceData() && actual.HasReferenceData()) {
        const PbiRawReferenceData& e2 = expected.ReferenceData();
        const PbiRawReferenceData& a2 = actual.ReferenceData();
        EXPECT_EQ(e2.entries_, a2.entries_);
    }

    // barcode data
    EXPECT_EQ(expected.HasBarcodeData(), actual.HasBarcodeData());
    if (expected.HasBarcodeData() && actual.HasBarcodeData()) {
        const PbiRawBarcodeData& e2 = expected.BarcodeData();
        const PbiRawBarcodeData& a2 = actual.BarcodeData();
        EXPECT_EQ(e2.bcForward_, a2.bcForward_);
        EXPECT_EQ(e2.bcReverse_, a2.bcReverse_);
        EXPECT_EQ(e2.bcQual_, a2.bcQual_);
    }
}

}  // namespace PacBioIndexTests

TEST(PacBioIndexTest, CreateFromExistingBam)
{
    // do this in temp directory, so we can ensure write access
    const std::string tempDir = PbbamTestsConfig::GeneratedData_Dir + "/";
    const std::string tempBamFn = tempDir + "aligned_copy.bam";
    const std::string tempPbiFn = tempBamFn + ".pbi";
    std::string cmd("cp ");
    cmd += PacBioIndexTests::test2BamFn;
    cmd += " ";
    cmd += tempBamFn;
    const auto cmdResult = system(cmd.c_str());
    UNUSED(cmdResult);

    BamFile bamFile(tempBamFn);
    PbiFile::CreateFrom(bamFile);
    EXPECT_EQ(tempPbiFn, bamFile.PacBioIndexFilename());

    PbiRawData index(bamFile.PacBioIndexFilename());
    EXPECT_EQ(PbiFile::Version_3_0_1, index.Version());
    EXPECT_EQ(10, index.NumReads());
    EXPECT_TRUE(index.HasMappedData());

    const PbiRawData expectedIndex = PacBioIndexTests::Test2Bam_ExistingIndex();
    PacBioIndexTests::ExpectRawIndicesEqual(expectedIndex, index);

    // clean up temp file(s)
    remove(tempBamFn.c_str());
    remove(tempPbiFn.c_str());
}

::testing::AssertionResult CanRead(BamReader& reader, BamRecord& record, int i)
{
    if (reader.GetNext(record))
        return ::testing::AssertionSuccess() << "i: " << i;
    else
        return ::testing::AssertionFailure() << "i: " << i;
}

TEST(PacBioIndexTest, CreateOnTheFly)
{
    // do this in temp directory, so we can ensure write access
    const std::string tempDir = PbbamTestsConfig::GeneratedData_Dir + "/";
    const std::string tempBamFn = tempDir + "temp.bam";
    const std::string tempPbiFn = tempBamFn + ".pbi";

    // NOTE: new file differs in size than existing (different write parameters may yield different file sizes, even though content is same)
    const std::vector<int64_t> expectedNewOffsets = {33816576,   236126208, 391315456, 469106688,
                                                     537067520,  587792384, 867303424, 1182793728,
                                                     1449787392, 1582628864};
    std::vector<int64_t> observedOffsets;

    // create PBI on the fly from input BAM while we write to new file
    {
        BamFile bamFile(PacBioIndexTests::test2BamFn);
        BamHeader header = bamFile.Header();

        BamWriter writer(tempBamFn, header);  // default compression, default thread count
        PbiBuilder builder(tempPbiFn, header.Sequences().size());

        int64_t vOffset = 0;
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile) {
            writer.Write(record, &vOffset);
            builder.AddRecord(record, vOffset);
            observedOffsets.push_back(vOffset);
        }
    }

    EXPECT_EQ(expectedNewOffsets, observedOffsets);

    // sanity check on original file
    {
        const std::vector<int64_t> originalFileOffsets = {33816576, 33825163,  33831333, 33834264,
                                                          33836542, 33838065,  33849818, 33863499,
                                                          33874621, 1392836608};
        BamRecord r;
        BamReader reader(PacBioIndexTests::test2BamFn);
        for (size_t i = 0; i < originalFileOffsets.size(); ++i) {
            reader.VirtualSeek(originalFileOffsets.at(i));
            EXPECT_TRUE(CanRead(reader, r, i));
        }
    }

    // attempt to seek in our new file using both expected & observed offsets
    {
        BamRecord r;
        BamReader reader(tempBamFn);
        for (size_t i = 0; i < expectedNewOffsets.size(); ++i) {
            reader.VirtualSeek(expectedNewOffsets.at(i));
            EXPECT_TRUE(CanRead(reader, r, i));
        }
        for (size_t i = 0; i < observedOffsets.size(); ++i) {
            reader.VirtualSeek(observedOffsets.at(i));
            EXPECT_TRUE(CanRead(reader, r, i));
        }
    }

    // compare data in new PBI file, to expected data
    const PbiRawData expectedIndex = PacBioIndexTests::Test2Bam_NewIndex();
    const PbiRawData fromBuilt = PbiRawData(tempPbiFn);
    PacBioIndexTests::ExpectRawIndicesEqual(expectedIndex, fromBuilt);

    // straight diff of newly-generated PBI file to existing PBI
    // TODO: Come back to this once pbindexump is in place.
    //       We can't exactly do this since file offsets may differ between 2 BAMs of differing compression levels.
    //       Should add some sort of BAM checksum based on contents, not just size, for this reason.
    //    const string pbiDiffCmd = string("diff -q ") + PacBioIndexTests::test2BamFn + ".pbi " + tempPbiFn;
    //    EXPECT_EQ(0, system(pbiDiffCmd.c_str()));

    // clean up temp file(s)
    remove(tempBamFn.c_str());
    remove(tempPbiFn.c_str());
}

TEST(PacBioIndexTest, RawLoadFromPbiFile)
{
    const BamFile bamFile(PacBioIndexTests::test2BamFn);
    const std::string pbiFilename = bamFile.PacBioIndexFilename();
    const PbiRawData loadedIndex(pbiFilename);

    const PbiRawData expectedIndex = PacBioIndexTests::Test2Bam_ExistingIndex();
    PacBioIndexTests::ExpectRawIndicesEqual(expectedIndex, loadedIndex);
}

TEST(PacBioIndexTest, BasicAndBarodeSectionsOnly)
{
    // do this in temp directory, so we can ensure write access
    const std::string tempDir = PbbamTestsConfig::GeneratedData_Dir + "/";
    const std::string tempBamFn = tempDir + "phi29.bam";
    const std::string tempPbiFn = tempBamFn + ".pbi";
    std::string cmd("cp ");
    cmd += PacBioIndexTests::phi29BamFn;
    cmd += " ";
    cmd += tempDir;
    const auto cmdResult = system(cmd.c_str());
    UNUSED(cmdResult);

    BamFile bamFile(tempBamFn);
    PbiFile::CreateFrom(bamFile);
    EXPECT_EQ(tempPbiFn, bamFile.PacBioIndexFilename());

    PbiRawData index(bamFile.PacBioIndexFilename());
    EXPECT_EQ(PbiFile::Version_3_0_1, index.Version());
    EXPECT_EQ(120, index.NumReads());
    EXPECT_FALSE(index.HasMappedData());
    EXPECT_TRUE(index.HasBarcodeData());

    const std::vector<int16_t> expectedBcForward{
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    const std::vector<int16_t> expectedBcReverse{
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    const std::vector<int8_t> expectedBcQuality{
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    const PbiRawBarcodeData& barcodeData = index.BarcodeData();
    EXPECT_EQ(expectedBcForward, barcodeData.bcForward_);
    EXPECT_EQ(expectedBcReverse, barcodeData.bcReverse_);
    EXPECT_EQ(expectedBcQuality, barcodeData.bcQual_);

    // clean up temp file(s)
    remove(tempBamFn.c_str());
    remove(tempPbiFn.c_str());
}

TEST(PacBioIndexTest, ReferenceDataNotLoadedOnUnsortedBam)
{
    BamFile bamFile(PacBioIndexTests::test2BamFn);
    PbiRawData raw(bamFile.PacBioIndexFilename());
    EXPECT_TRUE(raw.HasReferenceData());
}

TEST(PacBioIndexTest, LookupLoadFromFileOk)
{
    const uint32_t expectedNumReads = 10;
    const std::vector<int64_t> expectedOffsets{33816576, 33825163, 33831333, 33834264, 33836542,
                                               33838065, 33849818, 33863499, 33874621, 1392836608};

    EXPECT_NO_THROW({
        BamFile bamFile(PacBioIndexTests::test2BamFn);
        PbiRawData index(bamFile.PacBioIndexFilename());
        EXPECT_EQ(expectedNumReads, index.NumReads());
        EXPECT_EQ(expectedOffsets, index.BasicData().fileOffset_);
    });
}

TEST(PacBioIndexTest, ThrowOnNonExistentPbiFile)
{
    EXPECT_THROW(PbiRawData("does_not_exist.pbi"), std::exception);
}

TEST(PacBioIndexTest, ThrowOnNonPbiFile)
{
    // completely wrong format
    EXPECT_THROW(
        {
            const auto fastaFn = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";
            PbiRawData idx{fastaFn};
        },
        std::runtime_error);

    // BGZF file, but not PBI
    EXPECT_THROW(
        {
            const auto bamFn = PbbamTestsConfig::Data_Dir + "/ex2.bam";
            PbiRawData idx{bamFn};
        },
        std::runtime_error);
}

TEST(PacBioIndexTest, AggregatePBI)
{

    DataSet ds;
    ExternalResources& resources = ds.ExternalResources();
    resources.Add(BamFile{PbbamTestsConfig::Data_Dir +
                          "/aligned.bam"});  // 4 reads, BASIC | MAPPED | REFERENCE
    resources.Add(BamFile{PbbamTestsConfig::Data_Dir +
                          "/polymerase/production.subreads.bam"});  // 8 reads, BASIC | BARCODE
    resources.Add(BamFile{PbbamTestsConfig::Data_Dir +
                          "/polymerase/production_hq.hqregion.bam"});  // 1 read,  BASIC only

    const PbiRawData index{ds};
    const PbiRawBasicData& mergedBasicData = index.BasicData();
    const PbiRawBarcodeData& mergedBarcodeData = index.BarcodeData();
    const PbiRawMappedData& mergedMappedData = index.MappedData();

    const uint32_t expectedTotal = 13;  // 4 + 8 + 1

    // 'meta' info
    EXPECT_EQ(expectedTotal, index.NumReads());
    EXPECT_EQ(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::BARCODE, index.FileSections());
    EXPECT_TRUE(index.HasBarcodeData());
    EXPECT_TRUE(index.HasMappedData());
    EXPECT_FALSE(index.HasReferenceData());

    // file numbers
    EXPECT_EQ(0, mergedBasicData.fileNumber_.at(0));
    EXPECT_EQ(0, mergedBasicData.fileNumber_.at(1));
    EXPECT_EQ(0, mergedBasicData.fileNumber_.at(2));
    EXPECT_EQ(0, mergedBasicData.fileNumber_.at(3));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(4));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(5));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(6));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(7));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(8));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(9));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(10));
    EXPECT_EQ(1, mergedBasicData.fileNumber_.at(11));
    EXPECT_EQ(2, mergedBasicData.fileNumber_.at(12));

    // basic data
    EXPECT_EQ(0, mergedBasicData.qStart_.at(0));  // file 1
    EXPECT_EQ(0, mergedBasicData.qStart_.at(1));
    EXPECT_EQ(2659, mergedBasicData.qStart_.at(4));  // file 2
    EXPECT_EQ(3116, mergedBasicData.qStart_.at(5));
    EXPECT_EQ(2659, mergedBasicData.qStart_.at(12));  // file 3

    EXPECT_EQ(21102592, mergedBasicData.fileOffset_.at(0));  // file 1
    EXPECT_EQ(21102883, mergedBasicData.fileOffset_.at(1));
    EXPECT_EQ(19857408, mergedBasicData.fileOffset_.at(4));  // file 2
    EXPECT_EQ(19860696, mergedBasicData.fileOffset_.at(5));
    EXPECT_EQ(20054016, mergedBasicData.fileOffset_.at(12));  // file 3

    // mapped data
    EXPECT_EQ(60, mergedMappedData.mapQV_.at(0));  // file 1
    EXPECT_EQ(60, mergedMappedData.mapQV_.at(1));
    EXPECT_EQ(255, mergedMappedData.mapQV_.at(4));  // file 2
    EXPECT_EQ(255, mergedMappedData.mapQV_.at(5));
    EXPECT_EQ(255, mergedMappedData.mapQV_.at(12));  // file 3

    // barcode data
    EXPECT_EQ(-1, mergedBarcodeData.bcForward_.at(0));  // file 1
    EXPECT_EQ(-1, mergedBarcodeData.bcForward_.at(1));
    EXPECT_EQ(92, mergedBarcodeData.bcForward_.at(4));  // file 2
    EXPECT_EQ(92, mergedBarcodeData.bcForward_.at(5));
    EXPECT_EQ(-1, mergedBarcodeData.bcForward_.at(12));  // file 3
}
