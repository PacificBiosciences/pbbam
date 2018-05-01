// Author: Derek Barnett

#include <cstddef>
#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/virtual/VirtualPolymeraseCompositeReader.h>
#include <pbbam/virtual/VirtualPolymeraseReader.h>
#include <pbbam/virtual/ZmwReadStitcher.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace ZmwReadStitcherTests {

static void Compare(const BamRecord& b1, const BamRecord& b2)
{
    EXPECT_TRUE(b1.HasDeletionQV());
    EXPECT_TRUE(b1.HasDeletionTag());
    EXPECT_TRUE(b1.HasInsertionQV());
    EXPECT_TRUE(b1.HasMergeQV());
    EXPECT_TRUE(b1.HasSubstitutionQV());
    EXPECT_TRUE(b1.HasSubstitutionTag());
    EXPECT_TRUE(b1.HasLabelQV());
    EXPECT_TRUE(b1.HasAltLabelQV());
    EXPECT_TRUE(b1.HasAltLabelTag());
    EXPECT_TRUE(b1.HasPkmean());
    EXPECT_TRUE(b1.HasPkmid());
    EXPECT_TRUE(b1.HasPulseCall());
    EXPECT_TRUE(b1.HasIPD());
    EXPECT_TRUE(b1.HasPulseWidth());
    EXPECT_TRUE(b1.HasPrePulseFrames());
    EXPECT_TRUE(b1.HasPulseCallWidth());
    EXPECT_TRUE(b1.HasPulseMergeQV());

    EXPECT_TRUE(b2.HasDeletionQV());
    EXPECT_TRUE(b2.HasDeletionTag());
    EXPECT_TRUE(b2.HasInsertionQV());
    EXPECT_TRUE(b2.HasMergeQV());
    EXPECT_TRUE(b2.HasSubstitutionQV());
    EXPECT_TRUE(b2.HasSubstitutionTag());
    EXPECT_TRUE(b2.HasLabelQV());
    EXPECT_TRUE(b2.HasAltLabelQV());
    EXPECT_TRUE(b2.HasAltLabelTag());
    EXPECT_TRUE(b2.HasPkmean());
    EXPECT_TRUE(b2.HasPkmid());
    EXPECT_TRUE(b2.HasPulseCall());
    EXPECT_TRUE(b2.HasIPD());
    EXPECT_TRUE(b2.HasPulseWidth());
    EXPECT_TRUE(b2.HasPrePulseFrames());
    EXPECT_TRUE(b2.HasPulseCallWidth());
    EXPECT_TRUE(b2.HasPulseMergeQV());

    EXPECT_EQ(b1.FullName(), b2.FullName());
    EXPECT_EQ(b1.HoleNumber(), b2.HoleNumber());
    EXPECT_EQ(b1.NumPasses(), b2.NumPasses());
    EXPECT_EQ(b1.Sequence(), b2.Sequence());
    EXPECT_EQ(b1.Qualities(), b2.Qualities());
    EXPECT_EQ(b1.DeletionQV(), b2.DeletionQV());
    EXPECT_EQ(b1.DeletionTag(), b2.DeletionTag());
    EXPECT_EQ(b1.InsertionQV(), b2.InsertionQV());
    EXPECT_EQ(b1.MergeQV(), b2.MergeQV());
    EXPECT_EQ(b1.SubstitutionQV(), b2.SubstitutionQV());
    EXPECT_EQ(b1.SubstitutionTag(), b2.SubstitutionTag());
    EXPECT_EQ(b1.LabelQV(), b2.LabelQV());
    EXPECT_EQ(b1.AltLabelQV(), b2.AltLabelQV());
    EXPECT_EQ(b1.AltLabelTag(), b2.AltLabelTag());
    EXPECT_EQ(b1.Pkmean(), b2.Pkmean());
    EXPECT_EQ(b1.Pkmid(), b2.Pkmid());
    EXPECT_EQ(b1.PulseCall(), b2.PulseCall());
    EXPECT_EQ(b1.IPD(), b2.IPD());
    EXPECT_EQ(b1.PulseWidth(), b2.PulseWidth());
    EXPECT_EQ(b1.PrePulseFrames(), b2.PrePulseFrames());
    EXPECT_EQ(b1.PulseCallWidth(), b2.PulseCallWidth());
    EXPECT_EQ(b1.ReadGroup(), b2.ReadGroup());
    EXPECT_EQ(b1.PulseMergeQV(), b2.PulseMergeQV());
}

static size_t NumVirtualRecords(const std::string& primaryBamFn, const std::string& scrapsBamFn)
{
    ZmwReadStitcher stitcher(primaryBamFn, scrapsBamFn);
    size_t count = 0;
    while (stitcher.HasNext()) {
        const auto record = stitcher.Next();
        //        ()record;
        ++count;
    }
    return count;
}

}  // namespace ZmwReadStitcherTests

TEST(ZmwReadStitching, FromBams_NoFilter)
{
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/internal.scraps.bam");
    size_t count = 0;
    while (stitcher.HasNext()) {
        const auto record = stitcher.Next();
        //        ()record;
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(ZmwReadStitching, FromBams_Filtered)
{
    PbiFilter filter{PbiZmwFilter{100000}};  // setup to match DataSet w/ filter
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/internal.scraps.bam",
                             filter);
    size_t count = 0;
    while (stitcher.HasNext()) {
        const auto record = stitcher.Next();
        EXPECT_EQ(100000, record.HoleNumber());
        ++count;
    }
    EXPECT_EQ(1, count);
}

TEST(ZmwReadStitching, FromDataSet_NoFilter)
{
    // dataset contains these resources (subreads/scraps + hqregion/scraps BAMs)
    const std::string primaryFn1 =
        PbbamTestsConfig::Data_Dir + "/polymerase/production.subreads.bam";
    const std::string scrapsFn1 = PbbamTestsConfig::Data_Dir + "/polymerase/production.scraps.bam";
    const std::string primaryFn2 =
        PbbamTestsConfig::Data_Dir + "/polymerase/production_hq.hqregion.bam";
    const std::string scrapsFn2 =
        PbbamTestsConfig::Data_Dir + "/polymerase/production_hq.scraps.bam";
    const size_t numExpectedRecords =
        ZmwReadStitcherTests::NumVirtualRecords(primaryFn1, scrapsFn1) +
        ZmwReadStitcherTests::NumVirtualRecords(primaryFn2, scrapsFn2);

    const std::string datasetFn =
        PbbamTestsConfig::Data_Dir + "/polymerase/multiple_resources.subread.dataset.xml";

    DataSet ds{datasetFn};
    ZmwReadStitcher stitcher{ds};
    size_t numObservedRecords = 0;
    while (stitcher.HasNext()) {
        const auto record = stitcher.Next();
        //        ()record;
        ++numObservedRecords;
    }
    EXPECT_EQ(numExpectedRecords, numObservedRecords);
}

TEST(ZmwReadStitching, FromDataSet_Filtered)
{
    // dataset contains these resources (subreads/scraps + hqregion/scraps BAMs)
    const std::string primaryFn1 =
        PbbamTestsConfig::Data_Dir + "/polymerase/production.subreads.bam";
    const std::string scrapsFn1 = PbbamTestsConfig::Data_Dir + "/polymerase/production.scraps.bam";
    const std::string primaryFn2 = PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam";
    const std::string scrapsFn2 = PbbamTestsConfig::Data_Dir + "/polymerase/internal.scraps.bam";
    const std::string primaryFn3 =
        PbbamTestsConfig::Data_Dir + "/polymerase/production_hq.hqregion.bam";
    const std::string scrapsFn3 =
        PbbamTestsConfig::Data_Dir + "/polymerase/production_hq.scraps.bam";
    const size_t totalRecords = ZmwReadStitcherTests::NumVirtualRecords(primaryFn1, scrapsFn1) +
                                ZmwReadStitcherTests::NumVirtualRecords(primaryFn2, scrapsFn2) +
                                ZmwReadStitcherTests::NumVirtualRecords(primaryFn3, scrapsFn3);
    EXPECT_EQ(5, totalRecords);

    // our filter will remove the 2 "production" BAM pairs
    // using a ZMW filter that only the "internal" pair should pass
    const std::string datasetFn =
        PbbamTestsConfig::Data_Dir + "/polymerase/filtered_resources.subread.dataset.xml";

    DataSet ds{datasetFn};
    ZmwReadStitcher stitcher{ds};
    size_t numObservedRecords = 0;
    while (stitcher.HasNext()) {
        const auto record = stitcher.Next();
        //        ()record;
        ++numObservedRecords;
    }
    EXPECT_EQ(1, numObservedRecords);
}

TEST(ZmwReadStitching, FromDataSet_EmptyDataSet)
{
    ZmwReadStitcher stitcher{DataSet{}};
    EXPECT_FALSE(stitcher.HasNext());
}

TEST(ZmwReadStitching, EmptyScrapsFile)
{
    const std::string primaryBamFn =
        PbbamTestsConfig::Data_Dir + "/polymerase/scrapless.subreads.bam";
    const std::string scrapsBamFn = PbbamTestsConfig::Data_Dir + "/polymerase/scrapless.scraps.bam";

    const BamFile primaryBam(primaryBamFn);
    const BamFile scrapsBam(scrapsBamFn);
    const PbiRawData primaryIdx(primaryBam.PacBioIndexFilename());
    const PbiRawData scrapsIdx(scrapsBam.PacBioIndexFilename());
    EXPECT_EQ(3, primaryIdx.NumReads());
    EXPECT_EQ(0, scrapsIdx.NumReads());

    int count = 0;
    ZmwReadStitcher stitcher(primaryBamFn, scrapsBamFn);
    while (stitcher.HasNext()) {
        auto record = stitcher.Next();
        //        ()record;
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(ZmwReadStitching, VirtualRegions)
{
    // Create virtual polymerase read
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/internal.scraps.bam");
    auto virtualRecord = stitcher.Next();

    auto regionMap = virtualRecord.VirtualRegionsMap();
    auto adapter = virtualRecord.VirtualRegionsTable(VirtualRegionType::ADAPTER);

    // Compare different accessors to same source
    EXPECT_EQ(regionMap[VirtualRegionType::ADAPTER], adapter);

    // Compare to truth
    EXPECT_EQ(3047, adapter[0].beginPos);
    EXPECT_EQ(3095, adapter[0].endPos);
    EXPECT_EQ(3650, adapter[1].beginPos);
    EXPECT_EQ(3700, adapter[1].endPos);
    EXPECT_EQ(4289, adapter[2].beginPos);
    EXPECT_EQ(4335, adapter[2].endPos);
    EXPECT_EQ(4888, adapter[3].beginPos);
    EXPECT_EQ(4939, adapter[3].endPos);
    EXPECT_EQ(5498, adapter[4].beginPos);
    EXPECT_EQ(5546, adapter[4].endPos);
    EXPECT_EQ(6116, adapter[5].beginPos);
    EXPECT_EQ(6173, adapter[5].endPos);
    EXPECT_EQ(6740, adapter[6].beginPos);
    EXPECT_EQ(6790, adapter[6].endPos);

    auto barcode = virtualRecord.VirtualRegionsTable(VirtualRegionType::BARCODE);
    EXPECT_EQ(regionMap[VirtualRegionType::BARCODE], barcode);
    EXPECT_EQ(3025, barcode[0].beginPos);
    EXPECT_EQ(3047, barcode[0].endPos);
    EXPECT_EQ(3095, barcode[1].beginPos);
    EXPECT_EQ(3116, barcode[1].endPos);
    EXPECT_EQ(3628, barcode[2].beginPos);
    EXPECT_EQ(3650, barcode[2].endPos);
    EXPECT_EQ(3700, barcode[3].beginPos);
    EXPECT_EQ(3722, barcode[3].endPos);
    EXPECT_EQ(4267, barcode[4].beginPos);
    EXPECT_EQ(4289, barcode[4].endPos);
    EXPECT_EQ(4335, barcode[5].beginPos);
    EXPECT_EQ(4356, barcode[5].endPos);
    EXPECT_EQ(4864, barcode[6].beginPos);
    EXPECT_EQ(4888, barcode[6].endPos);
    EXPECT_EQ(4939, barcode[7].beginPos);
    EXPECT_EQ(4960, barcode[7].endPos);
    EXPECT_EQ(5477, barcode[8].beginPos);
    EXPECT_EQ(5498, barcode[8].endPos);
    EXPECT_EQ(5546, barcode[9].beginPos);
    EXPECT_EQ(5571, barcode[9].endPos);
    EXPECT_EQ(6087, barcode[10].beginPos);
    EXPECT_EQ(6116, barcode[10].endPos);
    EXPECT_EQ(6173, barcode[11].beginPos);
    EXPECT_EQ(6199, barcode[11].endPos);
    EXPECT_EQ(6719, barcode[12].beginPos);
    EXPECT_EQ(6740, barcode[12].endPos);
    EXPECT_EQ(6790, barcode[13].beginPos);
    EXPECT_EQ(6812, barcode[13].endPos);

    auto lqregion = virtualRecord.VirtualRegionsTable(VirtualRegionType::LQREGION);
    EXPECT_EQ(regionMap[VirtualRegionType::LQREGION], lqregion);
    EXPECT_EQ(0, lqregion[0].beginPos);
    EXPECT_EQ(2659, lqregion[0].endPos);
    EXPECT_EQ(7034, lqregion[1].beginPos);
    EXPECT_EQ(7035, lqregion[1].endPos);

    auto hqregion = virtualRecord.VirtualRegionsTable(VirtualRegionType::HQREGION);
    EXPECT_EQ(regionMap[VirtualRegionType::HQREGION], hqregion);
    EXPECT_EQ(2659, hqregion[0].beginPos);
    EXPECT_EQ(7034, hqregion[0].endPos);
}

TEST(ZmwReadStitching, InternalSubreadsToOriginal)
{
    // Create virtual polymerase read
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/internal.scraps.bam");
    EXPECT_TRUE(stitcher.HasNext());
    auto virtualRecord = stitcher.Next();

    // Read original polymerase read
    BamFile polyBam(PbbamTestsConfig::Data_Dir + "/polymerase/internal.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);
    auto begin = polyQuery.begin();
    auto end = polyQuery.end();
    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin;

    // check
    ZmwReadStitcherTests::Compare(polyRecord, virtualRecord);
}

TEST(ZmwReadStitching, InternalHQToOriginal)
{
    // Create virtual polymerase read
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/internal.hqregions.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/internal.lqregions.bam");
    EXPECT_TRUE(stitcher.HasNext());
    auto virtualRecord = stitcher.Next();

    // Read original polymerase read
    BamFile polyBam(PbbamTestsConfig::Data_Dir + "/polymerase/internal.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);
    auto begin = polyQuery.begin();
    auto end = polyQuery.end();
    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin;

    // check
    ZmwReadStitcherTests::Compare(polyRecord, virtualRecord);
}

TEST(ZmwReadStitching, ProductionSubreadsToOriginal)
{
    // Create virtual polymerase read
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/production.subreads.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/production.scraps.bam");

    EXPECT_TRUE(stitcher.HasNext());
    auto virtualRecord = stitcher.Next();
    EXPECT_FALSE(stitcher.HasNext());

    // Read original polymerase read
    BamFile polyBam(PbbamTestsConfig::Data_Dir + "/polymerase/production.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);

    auto begin = polyQuery.begin();
    auto end = polyQuery.end();
    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin;

    EXPECT_EQ(polyRecord.FullName(), virtualRecord.FullName());
    EXPECT_EQ(polyRecord.HoleNumber(), virtualRecord.HoleNumber());
    EXPECT_FLOAT_EQ(polyRecord.ReadAccuracy(), virtualRecord.ReadAccuracy());
    EXPECT_EQ(polyRecord.NumPasses(), virtualRecord.NumPasses());
    EXPECT_EQ(polyRecord.Sequence(), virtualRecord.Sequence());
    EXPECT_EQ(polyRecord.Qualities(), virtualRecord.Qualities());
    EXPECT_EQ(polyRecord.DeletionQV(), virtualRecord.DeletionQV());
    EXPECT_EQ(polyRecord.DeletionTag(), virtualRecord.DeletionTag());
    EXPECT_EQ(polyRecord.InsertionQV(), virtualRecord.InsertionQV());
    EXPECT_EQ(polyRecord.MergeQV(), virtualRecord.MergeQV());
    EXPECT_EQ(polyRecord.SubstitutionQV(), virtualRecord.SubstitutionQV());
    EXPECT_EQ(polyRecord.SubstitutionTag(), virtualRecord.SubstitutionTag());
    EXPECT_EQ(polyRecord.IPD(), virtualRecord.IPDV1Frames());
    EXPECT_EQ(polyRecord.ReadGroup(), virtualRecord.ReadGroup());
}

TEST(ZmwReadStitching, ProductionHQToOriginal)
{
    // Create virtual polymerase read
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/production_hq.hqregion.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/production_hq.scraps.bam");
    EXPECT_TRUE(stitcher.HasNext());
    auto virtualRecord = stitcher.Next();
    EXPECT_FALSE(stitcher.HasNext());

    // Read original polymerase read
    BamFile polyBam(PbbamTestsConfig::Data_Dir + "/polymerase/production.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);

    auto begin = polyQuery.begin();
    auto end = polyQuery.end();
    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin;

    EXPECT_FALSE(polyRecord.HasPulseCall());
    EXPECT_FALSE(virtualRecord.HasPulseCall());
    EXPECT_EQ(polyRecord.FullName(), virtualRecord.FullName());
    EXPECT_EQ(polyRecord.HoleNumber(), virtualRecord.HoleNumber());
    EXPECT_EQ(polyRecord.ReadAccuracy(), virtualRecord.ReadAccuracy());
    EXPECT_EQ(polyRecord.NumPasses(), virtualRecord.NumPasses());
    EXPECT_EQ(polyRecord.Sequence(), virtualRecord.Sequence());
    EXPECT_EQ(polyRecord.Qualities(), virtualRecord.Qualities());
    EXPECT_EQ(polyRecord.DeletionQV(), virtualRecord.DeletionQV());
    EXPECT_EQ(polyRecord.DeletionTag(), virtualRecord.DeletionTag());
    EXPECT_EQ(polyRecord.InsertionQV(), virtualRecord.InsertionQV());
    EXPECT_EQ(polyRecord.MergeQV(), virtualRecord.MergeQV());
    EXPECT_EQ(polyRecord.SubstitutionQV(), virtualRecord.SubstitutionQV());
    EXPECT_EQ(polyRecord.SubstitutionTag(), virtualRecord.SubstitutionTag());
    EXPECT_EQ(polyRecord.IPD(), virtualRecord.IPDV1Frames());
    EXPECT_EQ(polyRecord.ReadGroup(), virtualRecord.ReadGroup());

    EXPECT_TRUE(polyRecord.HasDeletionQV());
    EXPECT_TRUE(polyRecord.HasDeletionTag());
    EXPECT_TRUE(polyRecord.HasInsertionQV());
    EXPECT_TRUE(polyRecord.HasMergeQV());
    EXPECT_TRUE(polyRecord.HasSubstitutionQV());
    EXPECT_TRUE(polyRecord.HasSubstitutionTag());
    EXPECT_TRUE(polyRecord.HasIPD());
    EXPECT_FALSE(polyRecord.HasLabelQV());
    EXPECT_FALSE(polyRecord.HasAltLabelQV());
    EXPECT_FALSE(polyRecord.HasAltLabelTag());
    EXPECT_FALSE(polyRecord.HasPkmean());
    EXPECT_FALSE(polyRecord.HasPkmid());
    EXPECT_FALSE(polyRecord.HasPulseCall());
    EXPECT_FALSE(polyRecord.HasPulseWidth());
    EXPECT_FALSE(polyRecord.HasPrePulseFrames());
    EXPECT_FALSE(polyRecord.HasPulseCallWidth());

    EXPECT_TRUE(virtualRecord.HasDeletionQV());
    EXPECT_TRUE(virtualRecord.HasDeletionTag());
    EXPECT_TRUE(virtualRecord.HasInsertionQV());
    EXPECT_TRUE(virtualRecord.HasMergeQV());
    EXPECT_TRUE(virtualRecord.HasSubstitutionQV());
    EXPECT_TRUE(virtualRecord.HasSubstitutionTag());
    EXPECT_TRUE(virtualRecord.HasIPD());
    EXPECT_FALSE(virtualRecord.HasLabelQV());
    EXPECT_FALSE(virtualRecord.HasAltLabelQV());
    EXPECT_FALSE(virtualRecord.HasAltLabelTag());
    EXPECT_FALSE(virtualRecord.HasPkmean());
    EXPECT_FALSE(virtualRecord.HasPkmid());
    EXPECT_FALSE(virtualRecord.HasPulseCall());
    EXPECT_FALSE(virtualRecord.HasPulseWidth());
    EXPECT_FALSE(virtualRecord.HasPrePulseFrames());
    EXPECT_FALSE(virtualRecord.HasPulseCallWidth());
}

TEST(ZmwReadStitching, VirtualRecord_VirtualRegionsTable)
{
    ZmwReadStitcher stitcher(PbbamTestsConfig::Data_Dir + "/polymerase/production.subreads.bam",
                             PbbamTestsConfig::Data_Dir + "/polymerase/production.scraps.bam");
    EXPECT_TRUE(stitcher.HasNext());
    const auto virtualRecord = stitcher.Next();

    const auto subreads = virtualRecord.VirtualRegionsTable(VirtualRegionType::SUBREAD);
    const auto adapters = virtualRecord.VirtualRegionsTable(VirtualRegionType::ADAPTER);
    const auto hqRegions = virtualRecord.VirtualRegionsTable(VirtualRegionType::HQREGION);
    const auto lqRegions = virtualRecord.VirtualRegionsTable(VirtualRegionType::LQREGION);
    const auto barcodes = virtualRecord.VirtualRegionsTable(VirtualRegionType::BARCODE);
    const auto filtered = virtualRecord.VirtualRegionsTable(VirtualRegionType::FILTERED);

    EXPECT_FALSE(subreads.empty());
    EXPECT_FALSE(adapters.empty());
    EXPECT_FALSE(hqRegions.empty());
    EXPECT_FALSE(lqRegions.empty());
    EXPECT_FALSE(barcodes.empty());
    EXPECT_TRUE(filtered.empty());  // this type not present in this data
}

TEST(ZmwReadStitching, LegacyTypedefsOk)
{
    {
        VirtualPolymeraseReader reader(
            PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam",
            PbbamTestsConfig::Data_Dir + "/polymerase/internal.scraps.bam");
        size_t count = 0;
        while (reader.HasNext()) {
            const auto record = reader.Next();
            //            ()record;
            ++count;
        }
        EXPECT_EQ(3, count);
    }

    {
        VirtualPolymeraseCompositeReader reader{DataSet{}};
        EXPECT_FALSE(reader.HasNext());
    }
}
