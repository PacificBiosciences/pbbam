// Copyright (c) 2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Armin Töpfer

#ifdef PBBAM_TESTING
#define private public
#endif

#include <iostream>
#include <map>
#include <string>

#include <gtest/gtest.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/Frames.h>
#include <pbbam/virtual/VirtualPolymeraseReader.h>
#include <pbbam/virtual/ZmwWhitelistVirtualReader.h>

#include "TestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace PacBio {
namespace BAM {
namespace tests {

void Compare(const BamRecord& b1, const BamRecord& b2)
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

    EXPECT_EQ(b1.FullName(),        b2.FullName());
    EXPECT_EQ(b1.HoleNumber(),      b2.HoleNumber());
    EXPECT_EQ(b1.NumPasses(),       b2.NumPasses());
    EXPECT_EQ(b1.Sequence(),        b2.Sequence());
    EXPECT_EQ(b1.Qualities(),       b2.Qualities());
    EXPECT_EQ(b1.DeletionQV(),      b2.DeletionQV());
    EXPECT_EQ(b1.DeletionTag(),     b2.DeletionTag());
    EXPECT_EQ(b1.InsertionQV(),     b2.InsertionQV());
    EXPECT_EQ(b1.MergeQV(),         b2.MergeQV());
    EXPECT_EQ(b1.SubstitutionQV(),  b2.SubstitutionQV());
    EXPECT_EQ(b1.SubstitutionTag(), b2.SubstitutionTag());
    EXPECT_EQ(b1.LabelQV(),         b2.LabelQV());
    EXPECT_EQ(b1.AltLabelQV(),      b2.AltLabelQV());
    EXPECT_EQ(b1.AltLabelTag(),     b2.AltLabelTag());
    EXPECT_EQ(b1.Pkmean(),          b2.Pkmean());
    EXPECT_EQ(b1.Pkmid(),           b2.Pkmid());
    EXPECT_EQ(b1.PulseCall(),       b2.PulseCall());
    EXPECT_EQ(b1.IPD(),             b2.IPD());
    EXPECT_EQ(b1.PulseWidth(),      b2.PulseWidth());
    EXPECT_EQ(b1.PrePulseFrames(),  b2.PrePulseFrames());
    EXPECT_EQ(b1.PulseCallWidth(),  b2.PulseCallWidth());
    EXPECT_EQ(b1.ReadGroup(),       b2.ReadGroup());
    EXPECT_EQ(b1.PulseMergeQV(),    b2.PulseMergeQV());
}

} // namespace tests
} // namespace BAM
} // namespace PacBio

TEST(VirtualPolymeraseReader, InternalSubreadsToOriginal)
{
	// Create virtual polymerase read
    VirtualPolymeraseReader vpr(tests::Data_Dir + "/polymerase/internal.subreads.bam",
                                tests::Data_Dir + "/polymerase/internal.scraps.bam");
    EXPECT_TRUE(vpr.HasNext());
    auto virtualRecord = vpr.Next();
    EXPECT_FALSE(vpr.HasNext());

    // Read original polymerase read
    BamFile polyBam(tests::Data_Dir + "/polymerase/internal.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);

    auto begin = polyQuery.begin();
    auto end = polyQuery.begin();

    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin++;
	EXPECT_TRUE(begin == end); 

    tests::Compare(polyRecord, virtualRecord);
}

TEST(VirtualPolymeraseReader, InternalHQToOriginal)
{
	// Create virtual polymerase read
    VirtualPolymeraseReader vpr(tests::Data_Dir + "/polymerase/internal.hqregions.bam",
                                tests::Data_Dir + "/polymerase/internal.lqregions.bam");
    EXPECT_TRUE(vpr.HasNext());
    auto virtualRecord = vpr.Next();
    EXPECT_FALSE(vpr.HasNext());

    // Read original polymerase read
    BamFile polyBam(tests::Data_Dir + "/polymerase/internal.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);

    auto begin = polyQuery.begin();
    auto end = polyQuery.begin();

    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin++;
	EXPECT_TRUE(begin == end);    

    tests::Compare(polyRecord, virtualRecord);
}

TEST(VirtualPolymeraseReader, VirtualRegions)
{
	// Create virtual polymerase read
    VirtualPolymeraseReader vpr(tests::Data_Dir + "/polymerase/internal.subreads.bam",
                                tests::Data_Dir + "/polymerase/internal.scraps.bam");
    auto virtualRecord = vpr.Next();

    auto regionMap = virtualRecord.VirtualRegionsMap();
    auto adapter = virtualRecord.VirtualRegionsTable(VirtualRegionType::ADAPTER);

    // Compare different accessors to same source
    EXPECT_EQ(regionMap[VirtualRegionType::ADAPTER], adapter);

    // Compare to truth
    EXPECT_EQ(3047,adapter[0].beginPos);
    EXPECT_EQ(3095,adapter[0].endPos);
    EXPECT_EQ(3650,adapter[1].beginPos);
    EXPECT_EQ(3700,adapter[1].endPos);
    EXPECT_EQ(4289,adapter[2].beginPos);
    EXPECT_EQ(4335,adapter[2].endPos);
    EXPECT_EQ(4888,adapter[3].beginPos);
    EXPECT_EQ(4939,adapter[3].endPos);
    EXPECT_EQ(5498,adapter[4].beginPos);
    EXPECT_EQ(5546,adapter[4].endPos);
    EXPECT_EQ(6116,adapter[5].beginPos);
    EXPECT_EQ(6173,adapter[5].endPos);
    EXPECT_EQ(6740,adapter[6].beginPos);
    EXPECT_EQ(6790,adapter[6].endPos);


    auto barcode = virtualRecord.VirtualRegionsTable(VirtualRegionType::BARCODE);
    EXPECT_EQ(regionMap[VirtualRegionType::BARCODE], barcode);

    EXPECT_EQ(3025,barcode[0].beginPos);
    EXPECT_EQ(3047,barcode[0].endPos);
    EXPECT_EQ(3095,barcode[1].beginPos);
    EXPECT_EQ(3116,barcode[1].endPos);
    EXPECT_EQ(3628,barcode[2].beginPos);
    EXPECT_EQ(3650,barcode[2].endPos);
    EXPECT_EQ(3700,barcode[3].beginPos);
    EXPECT_EQ(3722,barcode[3].endPos);
    EXPECT_EQ(4267,barcode[4].beginPos);
    EXPECT_EQ(4289,barcode[4].endPos);
    EXPECT_EQ(4335,barcode[5].beginPos);
    EXPECT_EQ(4356,barcode[5].endPos);
    EXPECT_EQ(4864,barcode[6].beginPos);
    EXPECT_EQ(4888,barcode[6].endPos);
    EXPECT_EQ(4939,barcode[7].beginPos);
    EXPECT_EQ(4960,barcode[7].endPos);
    EXPECT_EQ(5477,barcode[8].beginPos);
    EXPECT_EQ(5498,barcode[8].endPos);
    EXPECT_EQ(5546,barcode[9].beginPos);
    EXPECT_EQ(5571,barcode[9].endPos);
    EXPECT_EQ(6087,barcode[10].beginPos);
    EXPECT_EQ(6116,barcode[10].endPos);
    EXPECT_EQ(6173,barcode[11].beginPos);
    EXPECT_EQ(6199,barcode[11].endPos);
    EXPECT_EQ(6719,barcode[12].beginPos);
    EXPECT_EQ(6740,barcode[12].endPos);
    EXPECT_EQ(6790,barcode[13].beginPos);
    EXPECT_EQ(6812,barcode[13].endPos);


    auto lqregion = virtualRecord.VirtualRegionsTable(VirtualRegionType::LQREGION);
    EXPECT_EQ(regionMap[VirtualRegionType::LQREGION], lqregion);

    EXPECT_EQ(0,lqregion[0].beginPos);
    EXPECT_EQ(2659,lqregion[0].endPos);
    EXPECT_EQ(7034,lqregion[1].beginPos);
    EXPECT_EQ(7035,lqregion[1].endPos);


    auto hqregion = virtualRecord.VirtualRegionsTable(VirtualRegionType::HQREGION);
    EXPECT_EQ(regionMap[VirtualRegionType::HQREGION], hqregion);

    EXPECT_EQ(2659,hqregion[0].beginPos);
    EXPECT_EQ(7034,hqregion[0].endPos);
}

TEST(VirtualPolymeraseReader, ProductionSubreadsToOriginal)
{
    // Create virtual polymerase read
    VirtualPolymeraseReader vpr(tests::Data_Dir + "/polymerase/production.subreads.bam",
                                tests::Data_Dir + "/polymerase/production.scraps.bam");

    EXPECT_TRUE(vpr.HasNext());
    auto virtualRecord = vpr.Next();
    EXPECT_FALSE(vpr.HasNext());

    // Read original polymerase read
    BamFile polyBam(tests::Data_Dir + "/polymerase/production.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);

    auto begin = polyQuery.begin();
    auto end = polyQuery.begin();

    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin++;
    EXPECT_TRUE(begin == end);

    EXPECT_EQ(polyRecord.FullName(),        virtualRecord.FullName());
    EXPECT_EQ(polyRecord.HoleNumber(),      virtualRecord.HoleNumber());
    EXPECT_FLOAT_EQ(polyRecord.ReadAccuracy(),    virtualRecord.ReadAccuracy());
    EXPECT_EQ(polyRecord.NumPasses(),       virtualRecord.NumPasses());
    EXPECT_EQ(polyRecord.Sequence(),        virtualRecord.Sequence());
    EXPECT_EQ(polyRecord.Qualities(),       virtualRecord.Qualities());
    EXPECT_EQ(polyRecord.DeletionQV(),      virtualRecord.DeletionQV());
    EXPECT_EQ(polyRecord.DeletionTag(),     virtualRecord.DeletionTag());
    EXPECT_EQ(polyRecord.InsertionQV(),     virtualRecord.InsertionQV());
    EXPECT_EQ(polyRecord.MergeQV(),         virtualRecord.MergeQV());
    EXPECT_EQ(polyRecord.SubstitutionQV(),  virtualRecord.SubstitutionQV());
    EXPECT_EQ(polyRecord.SubstitutionTag(), virtualRecord.SubstitutionTag());
    EXPECT_EQ(polyRecord.IPD(),             virtualRecord.IPDV1Frames());
    EXPECT_EQ(polyRecord.ReadGroup(),       virtualRecord.ReadGroup());
}

TEST(VirtualPolymeraseReader, ProductionHQToOriginal)
{
    // Create virtual polymerase read
    VirtualPolymeraseReader vpr(tests::Data_Dir + "/polymerase/production_hq.hqregion.bam",
                                tests::Data_Dir + "/polymerase/production_hq.scraps.bam");
    EXPECT_TRUE(vpr.HasNext());
    auto virtualRecord = vpr.Next();
    EXPECT_FALSE(vpr.HasNext());

    // Read original polymerase read
    BamFile polyBam(tests::Data_Dir + "/polymerase/production.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);

    auto begin = polyQuery.begin();
    auto end = polyQuery.begin();

    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin++;
    EXPECT_TRUE(begin == end); 

    EXPECT_FALSE(polyRecord.HasPulseCall());
    EXPECT_FALSE(virtualRecord.HasPulseCall());
    EXPECT_EQ(polyRecord.FullName(),        virtualRecord.FullName());
    EXPECT_EQ(polyRecord.HoleNumber(),      virtualRecord.HoleNumber());
    EXPECT_EQ(polyRecord.ReadAccuracy(),    virtualRecord.ReadAccuracy());
    EXPECT_EQ(polyRecord.NumPasses(),       virtualRecord.NumPasses());
    EXPECT_EQ(polyRecord.Sequence(),        virtualRecord.Sequence());
    EXPECT_EQ(polyRecord.Qualities(),       virtualRecord.Qualities());
    EXPECT_EQ(polyRecord.DeletionQV(),      virtualRecord.DeletionQV());
    EXPECT_EQ(polyRecord.DeletionTag(),     virtualRecord.DeletionTag());
    EXPECT_EQ(polyRecord.InsertionQV(),     virtualRecord.InsertionQV());
    EXPECT_EQ(polyRecord.MergeQV(),         virtualRecord.MergeQV());
    EXPECT_EQ(polyRecord.SubstitutionQV(),  virtualRecord.SubstitutionQV());
    EXPECT_EQ(polyRecord.SubstitutionTag(), virtualRecord.SubstitutionTag());
    EXPECT_EQ(polyRecord.IPD(),             virtualRecord.IPDV1Frames());
    EXPECT_EQ(polyRecord.ReadGroup(),       virtualRecord.ReadGroup());

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

TEST(ZmwWhitelistVirtualReader, SingleZmwOk)
{
    const std::vector<int32_t> whitelist = { 200000 };

    ZmwWhitelistVirtualReader reader(whitelist,
                                     tests::Data_Dir + "/polymerase/whitelist/internal.subreads.bam",
                                     tests::Data_Dir + "/polymerase/whitelist/internal.scraps.bam");

    // create virtual record
    EXPECT_TRUE(reader.HasNext());
    auto virtualRecord = reader.Next();
    EXPECT_FALSE(reader.HasNext());

    // fetch original polymerase read (2nd record)
    BamFile polyBam(tests::Data_Dir + "/polymerase/whitelist/internal.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);
    auto begin = polyQuery.begin();
    auto end = polyQuery.end();
    EXPECT_TRUE(begin != end);
    ++begin;
    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin++;

    EXPECT_EQ(200000, virtualRecord.HoleNumber());

    tests::Compare(polyRecord, virtualRecord);
}

TEST(ZmwWhitelistVirtualReader, MultiZmwsOk)
{
    const std::vector<int32_t> whitelist = { 100000, 300000 };

    ZmwWhitelistVirtualReader reader(whitelist,
                                     tests::Data_Dir + "/polymerase/whitelist/internal.subreads.bam",
                                     tests::Data_Dir + "/polymerase/whitelist/internal.scraps.bam");


    // create virtual records
    EXPECT_TRUE(reader.HasNext());
    auto virtualRecord1 = reader.Next();
    EXPECT_TRUE(reader.HasNext());
    auto virtualRecord2 = reader.Next();
    EXPECT_FALSE(reader.HasNext());

    // fetch original polymerase reads (2nd record)
    BamFile polyBam(tests::Data_Dir + "/polymerase/whitelist/internal.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);
    auto begin = polyQuery.begin();
    auto end = polyQuery.end();

    EXPECT_TRUE(begin != end);
    auto polyRecord1 = *begin++;
    EXPECT_TRUE(begin != end);
    ++begin;
    EXPECT_TRUE(begin != end);
    auto polyRecord2 = *begin++;
    EXPECT_TRUE(begin == end);

    EXPECT_EQ(100000, virtualRecord1.HoleNumber());
    EXPECT_EQ(300000, virtualRecord2.HoleNumber());

    tests::Compare(polyRecord1, virtualRecord1);
    tests::Compare(polyRecord2, virtualRecord2);
}

TEST(ZmwWhitelistVirtualReader, EmptyListOk)
{
    const std::vector<int32_t> whitelist = { };

    ZmwWhitelistVirtualReader reader(whitelist,
                                     tests::Data_Dir + "/polymerase/whitelist/internal.subreads.bam",
                                     tests::Data_Dir + "/polymerase/whitelist/internal.scraps.bam");
    EXPECT_FALSE(reader.HasNext());
    EXPECT_TRUE(reader.NextRaw().empty());
}

TEST(ZmwWhitelistVirtualReader, EmptyScrapsFileOk)
{
    const std::vector<int32_t> whitelist = { 10944689, 10944690 };
    const std::string primaryBamFn = tests::Data_Dir + "/polymerase/whitelist/scrapless.subreads.bam" ;
    const std::string scrapsBamFn  = tests::Data_Dir + "/polymerase/whitelist/scrapless.scraps.bam" ;

    int count = 0;
    ZmwWhitelistVirtualReader reader(whitelist, primaryBamFn, scrapsBamFn);
    while (reader.HasNext()) {
        auto record = reader.Next();
        (void)record;
        ++count;
    }
    EXPECT_EQ(2, count);

    const BamFile primaryBam(primaryBamFn);
    const BamFile scrapsBam(scrapsBamFn);
    const PbiRawData primaryIdx(primaryBam.PacBioIndexFilename());
    const PbiRawData scrapsIdx(scrapsBam.PacBioIndexFilename());
    EXPECT_EQ(3, primaryIdx.NumReads());
    EXPECT_EQ(0, scrapsIdx.NumReads());
}

TEST(ZmwWhitelistVirtualReader, UnknownZmwOk)
{
    const std::vector<int32_t> whitelist = { 42 }; // ZMW not in our files

    ZmwWhitelistVirtualReader reader(whitelist,
                                     tests::Data_Dir + "/polymerase/whitelist/internal.subreads.bam",
                                     tests::Data_Dir + "/polymerase/whitelist/internal.scraps.bam");
    EXPECT_FALSE(reader.HasNext());
    EXPECT_TRUE(reader.NextRaw().empty());
}

TEST(ZmwWhitelistVirtualReader, MixedKnownAndUnknownZmwsOk)
{
    const std::vector<int32_t> whitelist = { 42, 200000, 24 };

    ZmwWhitelistVirtualReader reader(whitelist,
                                     tests::Data_Dir + "/polymerase/whitelist/internal.subreads.bam",
                                     tests::Data_Dir + "/polymerase/whitelist/internal.scraps.bam");

    // everything below should behave exactly as 'SingleValueOk' test,
    // as the unknown ZMWs will have been removed during construction

    // create virtual record
    EXPECT_TRUE(reader.HasNext());
    auto virtualRecord = reader.Next();
    EXPECT_FALSE(reader.HasNext());

    // fetch original polymerase read (2nd record)
    BamFile polyBam(tests::Data_Dir + "/polymerase/whitelist/internal.polymerase.bam");
    EntireFileQuery polyQuery(polyBam);
    auto begin = polyQuery.begin();
    auto end = polyQuery.end();
    EXPECT_TRUE(begin != end);
    ++begin;
    EXPECT_TRUE(begin != end);
    auto polyRecord = *begin++;

    EXPECT_EQ(200000, virtualRecord.HoleNumber());

    tests::Compare(polyRecord, virtualRecord);
}

TEST(VirtualPolymeraseBamRecord, VirtualRegionsTableOk) 
{
    VirtualPolymeraseReader vpr(tests::Data_Dir + "/polymerase/production.subreads.bam",
                                tests::Data_Dir + "/polymerase/production.scraps.bam");
    EXPECT_TRUE(vpr.HasNext());
    const auto virtualRecord = vpr.Next();

    const auto subreads  = virtualRecord.VirtualRegionsTable(VirtualRegionType::SUBREAD);
    const auto adapters  = virtualRecord.VirtualRegionsTable(VirtualRegionType::ADAPTER);
    const auto hqRegions = virtualRecord.VirtualRegionsTable(VirtualRegionType::HQREGION);
    const auto lqRegions = virtualRecord.VirtualRegionsTable(VirtualRegionType::LQREGION);
    const auto barcodes  = virtualRecord.VirtualRegionsTable(VirtualRegionType::BARCODE);
    const auto filtered  = virtualRecord.VirtualRegionsTable(VirtualRegionType::FILTERED);

    EXPECT_FALSE(subreads.empty());
    EXPECT_FALSE(adapters.empty());
    EXPECT_FALSE(hqRegions.empty());
    EXPECT_FALSE(lqRegions.empty());
    EXPECT_FALSE(barcodes.empty());
    EXPECT_TRUE(filtered.empty());    // this annnotation type is not in data set
}
