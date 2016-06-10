// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#ifdef PBBAM_TESTING
#define private public
#endif

#include "TestData.h"
#include <gtest/gtest.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiBuilder.h>
#include <pbbam/PbiIndex.h>
#include <pbbam/PbiLookupData.h>
#include <pbbam/PbiRawData.h>
#include <string>
#include <cstdio>
#include <cstdlib>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string test2BamFn = tests::Data_Dir + "/aligned2.bam";
const string phi29BamFn = tests::Data_Dir + "/phi29.bam";

namespace PacBio {
namespace BAM {
namespace tests {

static
PbiRawData Test2Bam_CoreIndexData(void)

{
    PbiRawData rawData;
    rawData.Version(PbiFile::Version_3_0_1);
    rawData.FileSections(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::REFERENCE);
    rawData.NumReads(10);

    PbiRawBasicData& basicData = rawData.BasicData();
    basicData.rgId_       = { -1197849594, -1197849594, -1197849594, -1197849594, -1197849594, -1197849594, -1197849594, -1197849594, -1197849594, -1197849594 };
    basicData.qStart_     = {48,387,0,9936,10232,7468,5557,7285,426,7064};
    basicData.qEnd_       = {1132,1134,344,10187,10394,8906,7235,8657,1045,7421};
    basicData.holeNumber_ = {49050,32328,32328,6469,6469,30983,13473,13473,19915,30983};
    basicData.readQual_   = {0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6};
    basicData.ctxtFlag_   = {0,0,0,0,0,0,0,0,0,0};
    basicData.fileOffset_ = { 33816576, 33825163, 33831333, 33834264, 33836542, 33838065, 33849818, 33863499, 33874621, 1392836608 };

    PbiRawMappedData& mappedData = rawData.MappedData();
    mappedData.tId_       = {0,0,0,0,0,0,0,0,0,0};
    mappedData.tStart_    = {0,302,675,2170,2203,3572,4506,4507,4592,4669};
    mappedData.tEnd_      = {471,1019,1026,2397,2326,5015,6125,5850,5203,5011};
    mappedData.aStart_    = {653,395,1,9960,10271,7468,5574,7285,441,7075};
    mappedData.aEnd_      = {1129,1134,344,10185,10394,8906,7235,8647,1040,7418};
    mappedData.revStrand_ = {0,1,0,1,0,1,1,0,1,0};
    mappedData.nM_        = {460,704,339,216,118,1394,1581,1313,583,333};
    mappedData.nMM_       = {0,0,0,0,0,0,0,0,0,0};
    mappedData.mapQV_     = {254,254,254,254,254,254,254,254,254,254};

    PbiRawReferenceData& referenceData = rawData.ReferenceData();
    referenceData.entries_ = {
        PbiReferenceEntry{0,0,10},
        PbiReferenceEntry{4294967295,4294967295,4294967295}
    };

    return rawData;
}

// NOTE: We have 2 different sets of offsets because the copied, new file differs in size than the existing one.
//
//       Unsure which combination of write parameters were used on the original. Things like thread count,
//       compression level, etc. can effect compression ratio, BGZF block sizes, etc. even though the BAM record
//       content itself is equal. So we'll just track these index values separately, for now at least.
//
static
PbiRawData Test2Bam_ExistingIndex(void)
{
    PbiRawData index = Test2Bam_CoreIndexData();
    index.BasicData().fileOffset_ = { 33816576, 33825163, 33831333, 33834264, 33836542, 33838065, 33849818, 33863499, 33874621, 1392836608 };
    return index;
}

static
PbiRawData Test2Bam_NewIndex(void)
{
    PbiRawData index = Test2Bam_CoreIndexData();
    index.BasicData().fileOffset_ = { 33816576, 236126208, 391315456, 469106688, 537067520, 587792384, 867303424, 1182793728, 1449787392, 1582628864 };
    return index;
}

static
void ExpectRawIndicesEqual(const PbiRawData& expected, const PbiRawData& actual)
{
    // header data
    EXPECT_EQ(expected.Version(),      actual.Version());
    EXPECT_EQ(expected.FileSections(), actual.FileSections());
    EXPECT_EQ(expected.NumReads(),     actual.NumReads());

    // subread data
    const PbiRawBasicData& e = expected.BasicData();
    const PbiRawBasicData& a = actual.BasicData();
    EXPECT_EQ(e.rgId_,       a.rgId_);
    EXPECT_EQ(e.qStart_,     a.qStart_);
    EXPECT_EQ(e.qEnd_,       a.qEnd_);
    EXPECT_EQ(e.holeNumber_, a.holeNumber_);
    EXPECT_EQ(e.readQual_,   a.readQual_);
    EXPECT_EQ(e.ctxtFlag_,   a.ctxtFlag_);
    EXPECT_EQ(e.fileOffset_, a.fileOffset_);

    // mapped data
    EXPECT_EQ(expected.HasMappedData(), actual.HasMappedData());
    if (expected.HasMappedData() && actual.HasMappedData()) {
        const PbiRawMappedData& e = expected.MappedData();
        const PbiRawMappedData& a = actual.MappedData();
        EXPECT_EQ(e.tId_,       a.tId_);
        EXPECT_EQ(e.tStart_,    a.tStart_);
        EXPECT_EQ(e.tEnd_,      a.tEnd_);
        EXPECT_EQ(e.aStart_,    a.aStart_);
        EXPECT_EQ(e.aEnd_,      a.aEnd_);
        EXPECT_EQ(e.revStrand_, a.revStrand_);
        EXPECT_EQ(e.nM_,        a.nM_);
        EXPECT_EQ(e.nMM_,       a.nMM_);
        EXPECT_EQ(e.mapQV_,     a.mapQV_);
    }

    // reference data
    EXPECT_EQ(expected.HasReferenceData(), actual.HasReferenceData());
    if (expected.HasReferenceData() && actual.HasReferenceData()) {
        const PbiRawReferenceData& e = expected.ReferenceData();
        const PbiRawReferenceData& a = actual.ReferenceData();
        EXPECT_EQ(e.entries_, a.entries_);
    }

    // barcode data
    EXPECT_EQ(expected.HasBarcodeData(), actual.HasBarcodeData());
    if (expected.HasBarcodeData() && actual.HasBarcodeData()) {
        const PbiRawBarcodeData& e = expected.BarcodeData();
        const PbiRawBarcodeData& a = actual.BarcodeData();
        EXPECT_EQ(e.bcForward_,   a.bcForward_);
        EXPECT_EQ(e.bcReverse_,  a.bcReverse_);
        EXPECT_EQ(e.bcQual_,   a.bcQual_);
    }
}

static
bool BasicLookupsEqual(const BasicLookupData& lhs,
                         const BasicLookupData& rhs)
{
    return (lhs.rgId_ == rhs.rgId_ &&
            lhs.qStart_ == rhs.qStart_ &&
            lhs.qEnd_ == rhs.qEnd_ &&
            lhs.holeNumber_ == rhs.holeNumber_ &&
            lhs.readQual_ == rhs.readQual_ &&
            lhs.ctxtFlag_ == rhs.ctxtFlag_ &&
            lhs.fileOffset_ == rhs.fileOffset_);
}

static
bool MappedLookupsEqual(const MappedLookupData& lhs,
                        const MappedLookupData& rhs)
{
    return (lhs.tId_ == rhs.tId_ &&
            lhs.tStart_ == rhs.tStart_ &&
            lhs.tEnd_ == rhs.tEnd_ &&
            lhs.aStart_ == rhs.aStart_ &&
            lhs.aEnd_ == rhs.aEnd_ &&
            lhs.nM_ == rhs.nM_ &&
            lhs.nMM_ == rhs.nMM_ &&
            lhs.mapQV_ == rhs.mapQV_ &&
            lhs.forwardStrand_ == rhs.forwardStrand_ &&
            lhs.reverseStrand_ == rhs.reverseStrand_);
}

static
bool ReferenceLookupsEqual(const ReferenceLookupData& lhs,
                           const ReferenceLookupData& rhs)
{
    return lhs.references_ == rhs.references_;
}

static
bool BarcodeLookupsEqual(const BarcodeLookupData& lhs,
                         const BarcodeLookupData& rhs)
{
    return (lhs.bcForward_ == rhs.bcForward_ &&
            lhs.bcReverse_ == rhs.bcReverse_ &&
            lhs.bcQual_ == rhs.bcQual_);
}

static
bool PbiIndicesEqual(const PbiIndex& lhs, const PbiIndex& rhs)
{
    using namespace ::PacBio::BAM;
    const unique_ptr<internal::PbiIndexPrivate>& lhsImpl = lhs.d_;
    const unique_ptr<internal::PbiIndexPrivate>& rhsImpl = rhs.d_;
    if (lhsImpl == rhsImpl)
        return true;
    if (lhsImpl == nullptr || rhsImpl == nullptr)
        return false;

    // metadata compare
    if (lhsImpl->version_  != rhsImpl->version_  ||
        lhsImpl->sections_ != rhsImpl->sections_ ||
        lhsImpl->numReads_ != rhsImpl->numReads_)
    { return false; }

    // component compare
    if ( !BasicLookupsEqual(lhsImpl->basicData_,         rhsImpl->basicData_)   ||
         !MappedLookupsEqual(lhsImpl->mappedData_,       rhsImpl->mappedData_)    ||
         !ReferenceLookupsEqual(lhsImpl->referenceData_, rhsImpl->referenceData_) ||
         !BarcodeLookupsEqual(lhsImpl->barcodeData_,     rhsImpl->barcodeData_))
    { return false; }

    // if we get here, OK
    return true;
}

} // namespace tests
} // namespace BAM
} // namespace PacBio

TEST(PacBioIndexTest, CreateFromExistingBam)
{
    // do this in temp directory, so we can ensure write access
    const string tempDir    = "/tmp/";
    const string tempBamFn  = tempDir + "aligned2.bam";
    const string tempPbiFn  = tempBamFn + ".pbi";
    string cmd("cp ");
    cmd += test2BamFn;
    cmd += " ";
    cmd += tempDir;
    int cmdResult = system(cmd.c_str());
    (void)cmdResult;

    BamFile bamFile(tempBamFn);
    PbiFile::CreateFrom(bamFile);
    EXPECT_EQ(tempPbiFn, bamFile.PacBioIndexFilename());

    PbiRawData index(bamFile.PacBioIndexFilename());
    EXPECT_EQ(PbiFile::Version_3_0_1,  index.Version());
    EXPECT_EQ(10, index.NumReads());
    EXPECT_TRUE(index.HasMappedData());

    const PbiRawData& expectedIndex = tests::Test2Bam_ExistingIndex();
    tests::ExpectRawIndicesEqual(expectedIndex, index);

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
    const string tempDir    = "/tmp/";
    const string tempBamFn  = tempDir + "temp.bam";
    const string tempPbiFn  = tempBamFn + ".pbi";

    // NOTE: new file differs in size than existing (different write parameters may yield different file sizes, even though content is same)
    const vector<int64_t> expectedNewOffsets = { 33816576, 236126208, 391315456, 469106688, 537067520, 587792384, 867303424, 1182793728, 1449787392, 1582628864 };
    vector<int64_t> observedOffsets;

    // create PBI on the fly from input BAM while we write to new file
    {
        BamFile bamFile(test2BamFn);
        BamHeader header = bamFile.Header();

        BamWriter writer(tempBamFn, header); // default compression, default thread count
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
        const vector<int64_t> originalFileOffsets = { 33816576, 33825163, 33831333, 33834264, 33836542, 33838065, 33849818, 33863499, 33874621, 1392836608 };
        BamRecord r;
        BamReader reader(test2BamFn);
        for (int i = 0; i < originalFileOffsets.size(); ++i) {
            reader.VirtualSeek(originalFileOffsets.at(i));
            EXPECT_TRUE(CanRead(reader, r, i));
        }
    }

    // attempt to seek in our new file using both expected & observed offsets
    {
        BamRecord r;
        BamReader reader(tempBamFn);
        for (int i = 0; i < expectedNewOffsets.size(); ++i) {
            reader.VirtualSeek(expectedNewOffsets.at(i));
            EXPECT_TRUE(CanRead(reader, r, i));
        }
        for (int i = 0; i < observedOffsets.size(); ++i) {
            reader.VirtualSeek(observedOffsets.at(i));
            EXPECT_TRUE(CanRead(reader, r, i));
        }
    }

    // compare data in new PBI file, to expected data
    const PbiRawData& expectedIndex = tests::Test2Bam_NewIndex();
    const PbiRawData& fromBuilt = PbiRawData(tempPbiFn);
    tests::ExpectRawIndicesEqual(expectedIndex, fromBuilt);

    // straight diff of newly-generated PBI file to existing PBI
    // TODO: Come back to this once pbindexump is in place.
    //       We can't exactly do this since file offsets may differ between 2 BAMs of differing compression levels.
    //       Should add some sort of BAM checksum based on contents, not just size, for this reason.
//    const string pbiDiffCmd = string("diff -q ") + test2BamFn + ".pbi " + tempPbiFn;
//    EXPECT_EQ(0, system(pbiDiffCmd.c_str()));

    // clean up temp file(s)
    remove(tempBamFn.c_str());
    remove(tempPbiFn.c_str());
}

TEST(PacBioIndexTest, RawLoadFromPbiFile)
{
    const BamFile bamFile(test2BamFn);
    const string& pbiFilename = bamFile.PacBioIndexFilename();
    const PbiRawData loadedIndex(pbiFilename);

    const PbiRawData& expectedIndex = tests::Test2Bam_ExistingIndex();
    tests::ExpectRawIndicesEqual(expectedIndex, loadedIndex);
}

TEST(PacBioIndexTest, BasicAndBarodeSectionsOnly)
{
    // do this in temp directory, so we can ensure write access
    const string tempDir    = "/tmp/";
    const string tempBamFn  = tempDir + "phi29.bam";
    const string tempPbiFn  = tempBamFn + ".pbi";
    string cmd("cp ");
    cmd += phi29BamFn;
    cmd += " ";
    cmd += tempDir;
    int cmdResult = system(cmd.c_str());
    (void)cmdResult;

    BamFile bamFile(tempBamFn);
    PbiFile::CreateFrom(bamFile);
    EXPECT_EQ(tempPbiFn, bamFile.PacBioIndexFilename());

    PbiRawData index(bamFile.PacBioIndexFilename());
    EXPECT_EQ(PbiFile::Version_3_0_1,  index.Version());
    EXPECT_EQ(120, index.NumReads());
    EXPECT_FALSE(index.HasMappedData());
    EXPECT_TRUE(index.HasBarcodeData());

    const vector<int16_t> expectedBcForward = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
        2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
    const vector<int16_t> expectedBcReverse = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
        2,2,2,2,2,2,2,2,2,2,2,2,2,2};
    const vector<int8_t>  expectedBcQuality = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1};

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
    BamFile bamFile(test2BamFn);
    PbiRawData raw(bamFile.PacBioIndexFilename());
    EXPECT_TRUE(raw.HasReferenceData());
}

TEST(PacBioIndexTest, LookupLoadFromFileOk)
{
    BamFile bamFile(test2BamFn);
    EXPECT_NO_THROW(
    {
        PbiIndex index(bamFile.PacBioIndexFilename());
        EXPECT_EQ(10, index.NumReads());
        EXPECT_EQ(vector<int64_t>({ 33816576, 33825163, 33831333, 33834264, 33836542, 33838065, 33849818, 33863499, 33874621, 1392836608 }), index.BasicData().VirtualFileOffsets());
    });
}

TEST(PacBioIndexTest, ThrowOnNonExistentPbiFile)
{
    EXPECT_THROW(PbiRawData raw("does_not_exist.pbi"), std::exception);
    EXPECT_THROW(PbiIndex idx("does_not_exist.pbi"),   std::exception);
}

TEST(PacBioIndexTest, ThrowOnNonPbiFile)
{
    // completely wrong format
    const std::string fastaFn = tests::Data_Dir + "/lambdaNEB.fa";
    EXPECT_THROW(PbiRawData raw(fastaFn), std::exception);
    EXPECT_THROW(PbiIndex idx(fastaFn),   std::exception);

    // BGZF file, but not PBI
    const std::string& bamFn = tests::Data_Dir + "/ex2.bam";
    EXPECT_THROW(PbiRawData raw(bamFn), std::exception);
    EXPECT_THROW(PbiIndex idx(bamFn),   std::exception);
}

TEST(PacBioIndexTest, Copy_and_Move)
{
    const PbiIndex lookup(test2BamFn + ".pbi");

    const PbiIndex copyConstructed(lookup);

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    const PbiIndex moveConstructed(std::move(PbiIndex(test2BamFn + ".pbi")));
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    PbiIndex copyAssigned;
    copyAssigned = lookup;

    PbiIndex moveAssigned;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    moveAssigned = std::move(PbiIndex(test2BamFn + ".pbi"));
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, copyConstructed));
    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, moveConstructed));
    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, copyAssigned));
    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, moveAssigned));
}

TEST(PacBioIndexTest, OrderedLookup)
{
    using PacBio::BAM::IndexList;
    using PacBio::BAM::OrderedLookup;

    OrderedLookup<int>::container_type oRawData;
    oRawData[11] = { 0, 3, 4 };
    oRawData[20] = { 1 };
    oRawData[42] = { 2, 7, 8 };
    oRawData[10] = { 5 };
    oRawData[12] = { 6 };
    oRawData[99] = { 9 };

    OrderedLookup<int> oLookup(oRawData);

    // EQUAL
    EXPECT_EQ(IndexList({5}),       oLookup.LookupIndices(10, Compare::EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4}), oLookup.LookupIndices(11, Compare::EQUAL));
    EXPECT_EQ(IndexList({6}),       oLookup.LookupIndices(12, Compare::EQUAL));
    EXPECT_EQ(IndexList({1}),       oLookup.LookupIndices(20, Compare::EQUAL));
    EXPECT_EQ(IndexList({2, 7, 8}), oLookup.LookupIndices(42, Compare::EQUAL));
    EXPECT_EQ(IndexList({9}),       oLookup.LookupIndices(99, Compare::EQUAL));
    EXPECT_EQ(IndexList(),          oLookup.LookupIndices(66, Compare::EQUAL)); // does not exist

    // NOT_EQUAL
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 6, 7, 8, 9}),    oLookup.LookupIndices(10, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({1, 2, 5, 6, 7, 8, 9}),          oLookup.LookupIndices(11, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 7, 8, 9}),    oLookup.LookupIndices(12, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 2, 3, 4, 5, 6, 7, 8, 9}),    oLookup.LookupIndices(20, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 3, 4, 5, 6, 9}),          oLookup.LookupIndices(42, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8}),    oLookup.LookupIndices(99, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), oLookup.LookupIndices(66, Compare::NOT_EQUAL)); // does not exist

    // LESS_THAN
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), oLookup.LookupIndices(13, Compare::LESS_THAN));
    EXPECT_EQ(IndexList({0, 3, 4, 5}),    oLookup.LookupIndices(12, Compare::LESS_THAN));
    // do more checks

    // LESS_THAN_EQUAL
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), oLookup.LookupIndices(13, Compare::LESS_THAN_EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), oLookup.LookupIndices(12, Compare::LESS_THAN_EQUAL));
    // more checks?

    // GREATER_THAN
    EXPECT_EQ(IndexList({2,7,8,9}), oLookup.LookupIndices(41, Compare::GREATER_THAN));
    EXPECT_EQ(IndexList({9}),       oLookup.LookupIndices(42, Compare::GREATER_THAN));
    // more checks?

    // GREATER_THAN_EQUAL
    EXPECT_EQ(IndexList({2,7,8,9}), oLookup.LookupIndices(41, Compare::GREATER_THAN_EQUAL));
    EXPECT_EQ(IndexList({2,7,8,9}), oLookup.LookupIndices(42, Compare::GREATER_THAN_EQUAL));
    // more checks?
}

TEST(PacBioIndexTest, UnorderedLookup)
{
    using PacBio::BAM::IndexList;
    using PacBio::BAM::UnorderedLookup;

    UnorderedLookup<int>::container_type uRawData;
    uRawData[11] = { 0, 3, 4 };
    uRawData[20] = { 1 };
    uRawData[42] = { 2, 7, 8 };
    uRawData[10] = { 5 };
    uRawData[12] = { 6 };
    uRawData[99] = { 9 };

    UnorderedLookup<int> uLookup(uRawData);

    // EQUAL
    EXPECT_EQ(IndexList({5}),       uLookup.LookupIndices(10, Compare::EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4}), uLookup.LookupIndices(11, Compare::EQUAL));
    EXPECT_EQ(IndexList({6}),       uLookup.LookupIndices(12, Compare::EQUAL));
    EXPECT_EQ(IndexList({1}),       uLookup.LookupIndices(20, Compare::EQUAL));
    EXPECT_EQ(IndexList({2, 7, 8}), uLookup.LookupIndices(42, Compare::EQUAL));
    EXPECT_EQ(IndexList({9}),       uLookup.LookupIndices(99, Compare::EQUAL));
    EXPECT_EQ(IndexList(),          uLookup.LookupIndices(66, Compare::EQUAL)); // does not exist

    // NOT_EQUAL
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 6, 7, 8, 9}),    uLookup.LookupIndices(10, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({1, 2, 5, 6, 7, 8, 9}),          uLookup.LookupIndices(11, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 7, 8, 9}),    uLookup.LookupIndices(12, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 2, 3, 4, 5, 6, 7, 8, 9}),    uLookup.LookupIndices(20, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 3, 4, 5, 6, 9}),          uLookup.LookupIndices(42, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8}),    uLookup.LookupIndices(99, Compare::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), uLookup.LookupIndices(66, Compare::NOT_EQUAL)); // does not exist

    // LESS_THAN
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), uLookup.LookupIndices(13, Compare::LESS_THAN));
    EXPECT_EQ(IndexList({0, 3, 4, 5}),    uLookup.LookupIndices(12, Compare::LESS_THAN));
    // more checks?

    // LESS_THAN_EQUAL
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), uLookup.LookupIndices(13, Compare::LESS_THAN_EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), uLookup.LookupIndices(12, Compare::LESS_THAN_EQUAL));
    // more checks?

    // GREATER_THAN
    EXPECT_EQ(IndexList({2,7,8,9}), uLookup.LookupIndices(41, Compare::GREATER_THAN));
    EXPECT_EQ(IndexList({9}),       uLookup.LookupIndices(42, Compare::GREATER_THAN));
    // more checks?

    // GREATER_THAN_EQUAL
    EXPECT_EQ(uLookup.LookupIndices(41, Compare::GREATER_THAN_EQUAL), IndexList({2,7,8,9}));
    EXPECT_EQ(uLookup.LookupIndices(42, Compare::GREATER_THAN_EQUAL), IndexList({2,7,8,9}));
    // more checks?
}

TEST(PacBioIndexTest, MergeBlocks)
{
    using PacBio::BAM::IndexList;
    using PacBio::BAM::IndexResultBlock;
    using PacBio::BAM::IndexResultBlocks;
    using PacBio::BAM::mergedIndexBlocks;
    using PacBio::BAM::OrderedLookup;

    OrderedLookup<int>::container_type oRawData;
    oRawData[11] = { 0, 3, 4 };
    oRawData[20] = { 1 };
    oRawData[42] = { 2, 7, 8 };
    oRawData[10] = { 5 };
    oRawData[12] = { 6 };
    oRawData[99] = { 9 };

    OrderedLookup<int> oLookup(oRawData);

    // EQUAL
    auto mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(10, Compare::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(5, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(11, Compare::EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 1), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(3, 2), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(12, Compare::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(6, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(20, Compare::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(1, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(42, Compare::EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(2, 1), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(7, 2), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(99, Compare::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(9, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(66, Compare::EQUAL));
    EXPECT_TRUE(mergedBlocks.empty());

    // NOT_EQUAL
    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(10, Compare::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 5), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(6, 4), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(11, Compare::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(1, 2), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(5, 5), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(12, Compare::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 6), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(7, 3), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(20, Compare::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 1), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(2, 8), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(42, Compare::NOT_EQUAL));
    EXPECT_EQ(3, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 2), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(3, 4), mergedBlocks.at(1));
    EXPECT_EQ(IndexResultBlock(9, 1), mergedBlocks.at(2));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(99, Compare::NOT_EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 9), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(66, Compare::NOT_EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 10), mergedBlocks.at(0));
}

TEST(PacBioIndexTest, ApplyOffsetsToBlocks)
{
    using PacBio::BAM::BasicLookupData;
    using PacBio::BAM::IndexList;
    using PacBio::BAM::IndexResultBlock;
    using PacBio::BAM::IndexResultBlocks;
    using PacBio::BAM::mergedIndexBlocks;
    using PacBio::BAM::OrderedLookup;

    OrderedLookup<int>::container_type oRawData;
    oRawData[11] = { 0, 3, 4 };
    oRawData[20] = { 1 };
    oRawData[42] = { 2, 7, 8 };
    oRawData[10] = { 5 };
    oRawData[12] = { 6 };
    oRawData[99] = { 9 };

    OrderedLookup<int> oLookup(oRawData);
    auto mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(10, Compare::NOT_EQUAL));

    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 5), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(6, 4), mergedBlocks.at(1));

    BasicLookupData basicLookupData;
    basicLookupData.fileOffset_ = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    basicLookupData.ApplyOffsets(mergedBlocks);

    EXPECT_EQ(2,  mergedBlocks.size());
    EXPECT_EQ(0,  mergedBlocks.at(0).virtualOffset_);
    EXPECT_EQ(5,  mergedBlocks.at(0).numReads_);
    EXPECT_EQ(60, mergedBlocks.at(1).virtualOffset_);
    EXPECT_EQ(4,  mergedBlocks.at(1).numReads_);
}

TEST(PacBioIndexTest, LookupMulti)
{
    using PacBio::BAM::BasicLookupData;
    using PacBio::BAM::IndexList;
    using PacBio::BAM::IndexResultBlock;
    using PacBio::BAM::IndexResultBlocks;
    using PacBio::BAM::mergedIndexBlocks;
    using PacBio::BAM::UnorderedLookup;

    UnorderedLookup<int32_t>::container_type uRawData;
    uRawData[11] = { 0, 3, 4 };
    uRawData[20] = { 1 };
    uRawData[42] = { 2, 7, 8 };
    uRawData[10] = { 5 };
    uRawData[12] = { 6 };
    uRawData[99] = { 9 };

    BasicLookupData basicLookup;
    basicLookup.rgId_ = UnorderedLookup<int32_t>(uRawData);
    basicLookup.fileOffset_ = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };

    const std::vector<int32_t> whitelist = { 11, 42, 20 };
    const auto indices = basicLookup.IndicesMulti(BasicLookupData::RG_ID, whitelist);

    IndexResultBlocks mergedBlocks = mergedIndexBlocks(indices);
    basicLookup.ApplyOffsets(mergedBlocks);

    EXPECT_EQ(IndexList({0, 3, 4, 2, 7, 8, 1}), indices);
    EXPECT_EQ(2, mergedBlocks.size());

    const IndexResultBlock& block0 = mergedBlocks.at(0);
    EXPECT_EQ(0, block0.firstIndex_);
    EXPECT_EQ(5, block0.numReads_);
    EXPECT_EQ(0, block0.virtualOffset_);

    const IndexResultBlock& block1 = mergedBlocks.at(1);
    EXPECT_EQ(7,  block1.firstIndex_);
    EXPECT_EQ(2,  block1.numReads_);
    EXPECT_EQ(70, block1.virtualOffset_);
}

TEST(PacBioIndexTest, LookupAPI)
{
    const PbiIndex index(test2BamFn + ".pbi");
    const BasicLookupData& basicData = index.BasicData();
    const MappedLookupData& mappedData = index.MappedData();
    const BarcodeLookupData& barcodeData = index.BarcodeData();

    // rgId == x
    IndexResultBlocks rgResult = mergedIndexBlocks(basicData.Indices(BasicLookupData::RG_ID, -1197849594));
    basicData.ApplyOffsets(rgResult);
    EXPECT_EQ(1, rgResult.size());
    EXPECT_EQ(0, rgResult.at(0).firstIndex_);
    EXPECT_EQ(10, rgResult.at(0).numReads_);
    EXPECT_EQ(33816576, rgResult.at(0).virtualOffset_);

    // rg != x
    IndexResultBlocks notRgResult = mergedIndexBlocks(basicData.Indices(BasicLookupData::RG_ID,
                                                                        -1197849594,
                                                                        Compare::NOT_EQUAL));
    basicData.ApplyOffsets(notRgResult);
    EXPECT_TRUE(notRgResult.empty());

    // tEnd <= x
    IndexResultBlocks tEndLteResult = mergedIndexBlocks(mappedData.Indices(MappedLookupData::T_END,
                                                                            4500,
                                                                            Compare::LESS_THAN_EQUAL));
    basicData.ApplyOffsets(tEndLteResult);
    EXPECT_EQ(1, tEndLteResult.size());
    EXPECT_EQ(0, tEndLteResult.at(0).firstIndex_);
    EXPECT_EQ(5, tEndLteResult.at(0).numReads_);
    EXPECT_EQ(33816576, tEndLteResult.at(0).virtualOffset_);

    // tEnd >= x
    IndexResultBlocks tEndGteResult = mergedIndexBlocks(mappedData.Indices(MappedLookupData::T_START,
                                                                           4500,
                                                                           Compare::GREATER_THAN_EQUAL));
    basicData.ApplyOffsets(tEndGteResult);
    EXPECT_EQ(1, tEndGteResult.size());
    EXPECT_EQ(6, tEndGteResult.at(0).firstIndex_);
    EXPECT_EQ(4, tEndGteResult.at(0).numReads_);
    EXPECT_EQ(33849818, tEndGteResult.at(0).virtualOffset_);

    // strand query
    IndexResultBlocks forward = mergedIndexBlocks(mappedData.Indices(MappedLookupData::STRAND,
                                                                     Strand::FORWARD));
    basicData.ApplyOffsets(forward);
    EXPECT_EQ(5, forward.size());
    EXPECT_EQ(0, forward.at(0).firstIndex_);
    EXPECT_EQ(1, forward.at(0).numReads_);
    EXPECT_EQ(33816576, forward.at(0).virtualOffset_);

    EXPECT_EQ(2, forward.at(1).firstIndex_);
    EXPECT_EQ(1, forward.at(1).numReads_);
    EXPECT_EQ(33831333, forward.at(1).virtualOffset_);

    EXPECT_EQ(4, forward.at(2).firstIndex_);
    EXPECT_EQ(1, forward.at(2).numReads_);
    EXPECT_EQ(33836542, forward.at(2).virtualOffset_);

    EXPECT_EQ(7, forward.at(3).firstIndex_);
    EXPECT_EQ(1, forward.at(3).numReads_);
    EXPECT_EQ(33863499, forward.at(3).virtualOffset_);

    EXPECT_EQ(9, forward.at(4).firstIndex_);
    EXPECT_EQ(1, forward.at(4).numReads_);
    EXPECT_EQ(1392836608, forward.at(4).virtualOffset_);

    // 0,1,0,1,0,1,1,0,1,0
    IndexResultBlocks reverse = mergedIndexBlocks(mappedData.Indices(MappedLookupData::STRAND,
                                                                     Strand::REVERSE));
    basicData.ApplyOffsets(reverse);
    EXPECT_EQ(4, reverse.size());
    EXPECT_EQ(1, reverse.at(0).firstIndex_);
    EXPECT_EQ(1, reverse.at(0).numReads_);
    EXPECT_EQ(33825163, reverse.at(0).virtualOffset_);

    EXPECT_EQ(3, reverse.at(1).firstIndex_);
    EXPECT_EQ(1, reverse.at(1).numReads_);
    EXPECT_EQ(33834264, reverse.at(1).virtualOffset_);

    EXPECT_EQ(5, reverse.at(2).firstIndex_);
    EXPECT_EQ(2, reverse.at(2).numReads_);
    EXPECT_EQ(33838065, reverse.at(2).virtualOffset_);

    EXPECT_EQ(8, reverse.at(3).firstIndex_);
    EXPECT_EQ(1, reverse.at(3).numReads_);
    EXPECT_EQ(33874621, reverse.at(3).virtualOffset_);

    // query data field that is not in the PBI
    IndexResultBlocks missing = mergedIndexBlocks(barcodeData.Indices(BarcodeLookupData::BC_QUALITY,
                                                                      77,
                                                                      Compare::GREATER_THAN));
    basicData.ApplyOffsets(missing);
    EXPECT_TRUE(missing.empty());
}

TEST(PacBioIndexTest, LookupByZmw)
{
    BamFile f(tests::Data_Dir + "/dataset/bam_mapping.bam");
    f.EnsurePacBioIndexExists();

    const PbiIndex index(f.PacBioIndexFilename());
    const BasicLookupData& basicData = index.BasicData();

    IndexResultBlocks blocks =  mergedIndexBlocks(basicData.Indices(BasicLookupData::ZMW,
                                                                      20000,
                                                                      Compare::LESS_THAN));
    basicData.ApplyOffsets(blocks);
    EXPECT_EQ(14, blocks.size());

    //
    // we'll take a look at first 5 contiguous blocks of reads with ZMW < 20000
    //
    // skipped: { 49050, 32328, 32328 }
    // block0:  { 6469, 6469 }
    // skipped: { 30983 }
    // block1:  { 13473, 13473, 19915 }
    // skipped: { 30983 }
    // block2:  { 19915, 7247, 7247 }
    // skipped: { 38025 }
    // block3:  { 13473 }
    // skipped: { 36363, 36363, 31174, 31174, 38025, 50257, 50257 }
    // block4:  { 14743, 14743 }
    //

    const IndexResultBlock& block0 = blocks.at(0);
    EXPECT_EQ(3, block0.firstIndex_);
    EXPECT_EQ(2, block0.numReads_);
    EXPECT_EQ(32654529, block0.virtualOffset_);

    const IndexResultBlock& block1 = blocks.at(1);
    EXPECT_EQ(6, block1.firstIndex_);
    EXPECT_EQ(3, block1.numReads_);
    EXPECT_EQ(32669996, block1.virtualOffset_);

    const IndexResultBlock& block2 = blocks.at(2);
    EXPECT_EQ(10, block2.firstIndex_);
    EXPECT_EQ(3,  block2.numReads_);
    EXPECT_EQ(1388841957, block2.virtualOffset_);

    const IndexResultBlock& block3 = blocks.at(3);
    EXPECT_EQ(14, block3.firstIndex_);
    EXPECT_EQ(1,  block3.numReads_);
    EXPECT_EQ(1388864866, block3.virtualOffset_);

    const IndexResultBlock& block4 = blocks.at(4);
    EXPECT_EQ(22, block4.firstIndex_);
    EXPECT_EQ(2,  block4.numReads_);
    EXPECT_EQ(1388892121, block4.virtualOffset_);
}

TEST(PacBioIndexTest, LookupMultiZmw)
{
    BamFile f(tests::Data_Dir + "/dataset/bam_mapping.bam");
    f.EnsurePacBioIndexExists();

    const PbiIndex index(f.PacBioIndexFilename());
    const BasicLookupData& basicData = index.BasicData();

    const std::vector<int32_t> whitelist = { 13473, 38025 };
    IndexResultBlocks blocks = mergedIndexBlocks(basicData.IndicesMulti(BasicLookupData::ZMW, whitelist));
    basicData.ApplyOffsets(blocks);

    EXPECT_EQ(3, blocks.size());

    const IndexResultBlock& block0 = blocks.at(0);
    EXPECT_EQ(6, block0.firstIndex_);
    EXPECT_EQ(2, block0.numReads_);
    EXPECT_EQ(32669996, block0.virtualOffset_);

    const IndexResultBlock& block1 = blocks.at(1);
    EXPECT_EQ(13, block1.firstIndex_);
    EXPECT_EQ(2, block1.numReads_);
    EXPECT_EQ(1388851626, block1.virtualOffset_);

    const IndexResultBlock& block2 = blocks.at(2);
    EXPECT_EQ(19, block2.firstIndex_);
    EXPECT_EQ(1,  block2.numReads_);
    EXPECT_EQ(1388881468, block2.virtualOffset_);
}
