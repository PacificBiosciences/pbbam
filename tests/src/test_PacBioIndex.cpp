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
#include <pbbam/PbiIndex.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/internal/PbiIndex_p.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string test2BamFn = tests::Data_Dir + "/test_group_query/test2.bam";

namespace PacBio {
namespace BAM {

bool operator==(const PbiRawData& lhs, const PbiRawData& rhs)
{
    // header
    if (lhs.Version()      != rhs.Version())      return false;
    if (lhs.FileSections() != rhs.FileSections()) return false;
    if (lhs.NumReads()     != rhs.NumReads())     return false;

    // subread data
    const PbiRawSubreadData& lhsSubreads = lhs.SubreadData();
    const PbiRawSubreadData& rhsSubreads = rhs.SubreadData();
    if (lhsSubreads.rgId_       != rhsSubreads.rgId_)       return false;
    if (lhsSubreads.qStart_     != rhsSubreads.qStart_)     return false;
    if (lhsSubreads.qEnd_       != rhsSubreads.qEnd_)       return false;
    if (lhsSubreads.holeNumber_ != rhsSubreads.holeNumber_) return false;
    if (lhsSubreads.readQual_   != rhsSubreads.readQual_)   return false;
    if (lhsSubreads.fileOffset_ != rhsSubreads.fileOffset_) return false;

    // mapped data
    if (lhs.HasMappedData()) {
        if (!rhs.HasMappedData())
            return false;
        const PbiRawMappedData& lhsMapped = lhs.MappedData();
        const PbiRawMappedData& rhsMapped = rhs.MappedData();
        if (lhsMapped.tId_       != rhsMapped.tId_)       return false;
        if (lhsMapped.tStart_    != rhsMapped.tStart_)    return false;
        if (lhsMapped.tEnd_      != rhsMapped.tEnd_)      return false;
        if (lhsMapped.aStart_    != rhsMapped.aStart_)    return false;
        if (lhsMapped.aEnd_      != rhsMapped.aEnd_)      return false;
        if (lhsMapped.revStrand_ != rhsMapped.revStrand_) return false;
        if (lhsMapped.nM_        != rhsMapped.nM_)        return false;
        if (lhsMapped.nMM_       != rhsMapped.nMM_)       return false;
        if (lhsMapped.mapQV_     != rhsMapped.mapQV_)     return false;
    }

    // reference data
    if (lhs.HasReferenceData()) {
        if (!rhs.HasReferenceData())
            return false;
        const PbiRawReferenceData& lhsRefData = lhs.ReferenceData();
        const PbiRawReferenceData& rhsRefData = rhs.ReferenceData();
        if (lhsRefData.entries_ != rhsRefData.entries_ ) return false;
    }

    // barcode data
    if (lhs.HasBarcodeData()) {
        if (!rhs.HasBarcodeData())
            return false;
        const PbiRawBarcodeData& lhsBarcode = lhs.BarcodeData();
        const PbiRawBarcodeData& rhsBarcode = rhs.BarcodeData();
        if (lhsBarcode.bcLeft_   != rhsBarcode.bcLeft_)   return false;
        if (lhsBarcode.bcRight_  != rhsBarcode.bcRight_)  return false;
        if (lhsBarcode.bcQual_   != rhsBarcode.bcQual_)   return false;
        if (lhsBarcode.ctxtFlag_ != rhsBarcode.ctxtFlag_) return false;
    }

    // if we get here, OK
    return true;
}

namespace tests {

static
PbiRawData Test2Bam_RawIndex(void)
{
    PbiRawData index;
    index.Version(PbiFile::Version_3_0_0);
    index.FileSections(PbiFile::SUBREAD | PbiFile::MAPPED);
    index.NumReads(4);

    PbiRawSubreadData& subreadData = index.SubreadData();
    subreadData.rgId_       = { -1197849594, -1197849594, -1197849594, -1197849594 };
    subreadData.qStart_     = { 2114, 2579, 4101, 5615 };
    subreadData.qEnd_       = { 2531, 4055, 5571, 6237 };
    subreadData.holeNumber_ = { 14743, 14743, 14743, 14743 };
    subreadData.readQual_   = { 901, 901, 901, 901 };
    subreadData.fileOffset_ = { 35323904,  35327454, 35339466, 35351517 };

    PbiRawMappedData& mappedData = index.mappedData_;
    mappedData.tId_       = { 0, 0, 0, 0 };
    mappedData.tStart_    = { 9507, 8453, 8455, 9291 };
    mappedData.tEnd_      = { 9903, 9902, 9893, 9900 };
    mappedData.aStart_    = { 2130, 2130, 2130, 2130 };
    mappedData.aEnd_      = { 2531, 2531, 2531, 2531 };
    mappedData.revStrand_ = { 0, 1, 0, 1 };
    mappedData.nM_        = { 384, 1411, 1393, 598 };
    mappedData.nMM_       = { 0, 0, 0, 0 };
    mappedData.mapQV_     = { 254, 254, 254, 254 };

    // reference & barcode data are empty for this file
    return index;
}

static
bool SubreadLookupsEqual(const internal::SubreadLookupData& lhs,
                         const internal::SubreadLookupData& rhs)
{
    return (lhs.rgId_ == rhs.rgId_ &&
            lhs.qStart_ == rhs.qStart_ &&
            lhs.qEnd_ == rhs.qEnd_ &&
            lhs.holeNumber_ == rhs.holeNumber_ &&
            lhs.readQual_ == rhs.readQual_ &&
            lhs.fileOffset_ == rhs.fileOffset_);
}

static
bool MappedLookupsEqual(const internal::MappedLookupData& lhs,
                        const internal::MappedLookupData& rhs)
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
bool ReferenceLookupsEqual(const internal::ReferenceLookupData& lhs,
                           const internal::ReferenceLookupData& rhs)
{
    return lhs.references_ == rhs.references_;
}

static
bool BarcodeLookupsEqual(const internal::BarcodeLookupData& lhs,
                         const internal::BarcodeLookupData& rhs)
{
    return (lhs.bcLeft_ == rhs.bcLeft_ &&
            lhs.bcRight_ == rhs.bcRight_ &&
            lhs.bcQual_ == rhs.bcQual_ &&
            lhs.ctxtFlag_ == rhs.ctxtFlag_);
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
    if ( !SubreadLookupsEqual(lhsImpl->subreadData_,     rhsImpl->subreadData_)   ||
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

TEST(PacBioIndexTest, RawLoadFromFileOk)
{
    const BamFile bamFile(test2BamFn);
    const string& pbiFilename = bamFile.PacBioIndexFilename();
    const PbiRawData loadedIndex(pbiFilename);
    const PbiRawData& expectedIndex = tests::Test2Bam_RawIndex();
    EXPECT_TRUE(expectedIndex == loadedIndex);
}

TEST(PacBioIndexTest, BuildFromBamOk)
{
    BamFile bamFile(test2BamFn);
    PbiFile::CreateFrom(bamFile);

    PbiRawData index(bamFile.PacBioIndexFilename());
    EXPECT_EQ(PbiFile::Version_3_0_0,  index.Version());
    EXPECT_EQ(4, index.NumReads());
    EXPECT_TRUE(index.HasMappedData());

    const PbiRawData& expectedIndex = tests::Test2Bam_RawIndex();
    EXPECT_TRUE(expectedIndex == index);
}

TEST(PacBioIndexTest, ReferenceDataNotLoadedOnUnsortedBam)
{
    BamFile bamFile(test2BamFn);
    PbiFile::CreateFrom(bamFile);
    PbiRawData raw(bamFile.PacBioIndexFilename());
    EXPECT_FALSE(raw.HasReferenceData());
}

TEST(PacBioIndexTest, LookupLoadFromFileOk)
{
    BamFile bamFile(test2BamFn);
    EXPECT_NO_THROW(
    {
        PbiIndex index(bamFile.PacBioIndexFilename());
        EXPECT_EQ(4, index.NumReads());
        EXPECT_EQ(vector<int64_t>({ 35323904,  35327454, 35339466, 35351517 }), index.VirtualFileOffsets());
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
    const PbiIndex moveConstructed(std::move(PbiIndex(test2BamFn + ".pbi")));

    PbiIndex copyAssigned;
    copyAssigned = lookup;

    PbiIndex moveAssigned;
    moveAssigned = std::move(PbiIndex(test2BamFn + ".pbi"));

    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, copyConstructed));
    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, moveConstructed));
    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, copyAssigned));
    EXPECT_TRUE(tests::PbiIndicesEqual(lookup, moveAssigned));
}

TEST(PacBioIndexTest, OrderedLookup)
{
    using PacBio::BAM::CompareType;
    using PacBio::BAM::IndexList;
    using PacBio::BAM::internal::OrderedLookup;

    OrderedLookup<int>::ContainerType oRawData;
    oRawData[11] = { 0, 3, 4 };
    oRawData[20] = { 1 };
    oRawData[42] = { 2, 7, 8 };
    oRawData[10] = { 5 };
    oRawData[12] = { 6 };
    oRawData[99] = { 9 };

    OrderedLookup<int> oLookup(oRawData);

    // EQUAL
    EXPECT_EQ(IndexList({5}),       oLookup.LookupIndices(10, CompareType::EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4}), oLookup.LookupIndices(11, CompareType::EQUAL));
    EXPECT_EQ(IndexList({6}),       oLookup.LookupIndices(12, CompareType::EQUAL));
    EXPECT_EQ(IndexList({1}),       oLookup.LookupIndices(20, CompareType::EQUAL));
    EXPECT_EQ(IndexList({2, 7, 8}), oLookup.LookupIndices(42, CompareType::EQUAL));
    EXPECT_EQ(IndexList({9}),       oLookup.LookupIndices(99, CompareType::EQUAL));
    EXPECT_EQ(IndexList(),          oLookup.LookupIndices(66, CompareType::EQUAL)); // does not exist

    // NOT_EQUAL
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 6, 7, 8, 9}),    oLookup.LookupIndices(10, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({1, 2, 5, 6, 7, 8, 9}),          oLookup.LookupIndices(11, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 7, 8, 9}),    oLookup.LookupIndices(12, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 2, 3, 4, 5, 6, 7, 8, 9}),    oLookup.LookupIndices(20, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 3, 4, 5, 6, 9}),          oLookup.LookupIndices(42, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8}),    oLookup.LookupIndices(99, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), oLookup.LookupIndices(66, CompareType::NOT_EQUAL)); // does not exist

    // LESS_THAN
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), oLookup.LookupIndices(13, CompareType::LESS_THAN));
    EXPECT_EQ(IndexList({0, 3, 4, 5}),    oLookup.LookupIndices(12, CompareType::LESS_THAN));
    // do more checks

    // LESS_THAN_EQUAL
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), oLookup.LookupIndices(13, CompareType::LESS_THAN_EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), oLookup.LookupIndices(12, CompareType::LESS_THAN_EQUAL));
    // more checks?

    // GREATER_THAN
    EXPECT_EQ(IndexList({2,7,8,9}), oLookup.LookupIndices(41, CompareType::GREATER_THAN));
    EXPECT_EQ(IndexList({9}),       oLookup.LookupIndices(42, CompareType::GREATER_THAN));
    // more checks?

    // GREATER_THAN_EQUAL
    EXPECT_EQ(IndexList({2,7,8,9}), oLookup.LookupIndices(41, CompareType::GREATER_THAN_EQUAL));
    EXPECT_EQ(IndexList({2,7,8,9}), oLookup.LookupIndices(42, CompareType::GREATER_THAN_EQUAL));
    // more checks?
}

TEST(PacBioIndexTest, UnorderedLookup)
{
    using PacBio::BAM::CompareType;
    using PacBio::BAM::IndexList;
    using PacBio::BAM::internal::UnorderedLookup;

    UnorderedLookup<int>::ContainerType uRawData;
    uRawData[11] = { 0, 3, 4 };
    uRawData[20] = { 1 };
    uRawData[42] = { 2, 7, 8 };
    uRawData[10] = { 5 };
    uRawData[12] = { 6 };
    uRawData[99] = { 9 };

    UnorderedLookup<int> uLookup(uRawData);

    // EQUAL
    EXPECT_EQ(IndexList({5}),       uLookup.LookupIndices(10, CompareType::EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4}), uLookup.LookupIndices(11, CompareType::EQUAL));
    EXPECT_EQ(IndexList({6}),       uLookup.LookupIndices(12, CompareType::EQUAL));
    EXPECT_EQ(IndexList({1}),       uLookup.LookupIndices(20, CompareType::EQUAL));
    EXPECT_EQ(IndexList({2, 7, 8}), uLookup.LookupIndices(42, CompareType::EQUAL));
    EXPECT_EQ(IndexList({9}),       uLookup.LookupIndices(99, CompareType::EQUAL));
    EXPECT_EQ(IndexList(),          uLookup.LookupIndices(66, CompareType::EQUAL)); // does not exist

    // NOT_EQUAL
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 6, 7, 8, 9}),    uLookup.LookupIndices(10, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({1, 2, 5, 6, 7, 8, 9}),          uLookup.LookupIndices(11, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 7, 8, 9}),    uLookup.LookupIndices(12, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 2, 3, 4, 5, 6, 7, 8, 9}),    uLookup.LookupIndices(20, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 3, 4, 5, 6, 9}),          uLookup.LookupIndices(42, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8}),    uLookup.LookupIndices(99, CompareType::NOT_EQUAL));
    EXPECT_EQ(IndexList({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), uLookup.LookupIndices(66, CompareType::NOT_EQUAL)); // does not exist

    // LESS_THAN
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), uLookup.LookupIndices(13, CompareType::LESS_THAN));
    EXPECT_EQ(IndexList({0, 3, 4, 5}),    uLookup.LookupIndices(12, CompareType::LESS_THAN));
    // more checks?

    // LESS_THAN_EQUAL
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), uLookup.LookupIndices(13, CompareType::LESS_THAN_EQUAL));
    EXPECT_EQ(IndexList({0, 3, 4, 5, 6}), uLookup.LookupIndices(12, CompareType::LESS_THAN_EQUAL));
    // more checks?

    // GREATER_THAN
    EXPECT_EQ(IndexList({2,7,8,9}), uLookup.LookupIndices(41, CompareType::GREATER_THAN));
    EXPECT_EQ(IndexList({9}),       uLookup.LookupIndices(42, CompareType::GREATER_THAN));
    // more checks?

    // GREATER_THAN_EQUAL
    EXPECT_EQ(uLookup.LookupIndices(41, CompareType::GREATER_THAN_EQUAL), IndexList({2,7,8,9}));
    EXPECT_EQ(uLookup.LookupIndices(42, CompareType::GREATER_THAN_EQUAL), IndexList({2,7,8,9}));
    // more checks?
}

TEST(PacBioIndexTest, MergeBlocks)
{
    using PacBio::BAM::CompareType;
    using PacBio::BAM::IndexList;
    using PacBio::BAM::IndexResultBlock;
    using PacBio::BAM::IndexResultBlocks;
    using PacBio::BAM::internal::mergedIndexBlocks;
    using PacBio::BAM::internal::OrderedLookup;

    OrderedLookup<int>::ContainerType oRawData;
    oRawData[11] = { 0, 3, 4 };
    oRawData[20] = { 1 };
    oRawData[42] = { 2, 7, 8 };
    oRawData[10] = { 5 };
    oRawData[12] = { 6 };
    oRawData[99] = { 9 };

    OrderedLookup<int> oLookup(oRawData);

    // EQUAL
    auto mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(10, CompareType::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(5, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(11, CompareType::EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 1), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(3, 2), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(12, CompareType::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(6, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(20, CompareType::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(1, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(42, CompareType::EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(2, 1), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(7, 2), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(99, CompareType::EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(9, 1), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(66, CompareType::EQUAL));
    EXPECT_TRUE(mergedBlocks.empty());

    // NOT_EQUAL
    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(10, CompareType::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 5), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(6, 4), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(11, CompareType::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(1, 2), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(5, 5), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(12, CompareType::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 6), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(7, 3), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(20, CompareType::NOT_EQUAL));
    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 1), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(2, 8), mergedBlocks.at(1));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(42, CompareType::NOT_EQUAL));
    EXPECT_EQ(3, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 2), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(3, 4), mergedBlocks.at(1));
    EXPECT_EQ(IndexResultBlock(9, 1), mergedBlocks.at(2));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(99, CompareType::NOT_EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 9), mergedBlocks.at(0));

    mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(66, CompareType::NOT_EQUAL));
    EXPECT_EQ(1, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 10), mergedBlocks.at(0));
}

TEST(PacBioIndexTest, ApplyOffsetsToBlocks)
{
    using PacBio::BAM::CompareType;
    using PacBio::BAM::IndexList;
    using PacBio::BAM::IndexResultBlock;
    using PacBio::BAM::IndexResultBlocks;
    using PacBio::BAM::internal::mergedIndexBlocks;
    using PacBio::BAM::internal::OrderedLookup;
    using PacBio::BAM::internal::SubreadLookupData;

    OrderedLookup<int>::ContainerType oRawData;
    oRawData[11] = { 0, 3, 4 };
    oRawData[20] = { 1 };
    oRawData[42] = { 2, 7, 8 };
    oRawData[10] = { 5 };
    oRawData[12] = { 6 };
    oRawData[99] = { 9 };

    OrderedLookup<int> oLookup(oRawData);
    auto mergedBlocks = mergedIndexBlocks(oLookup.LookupIndices(10, CompareType::NOT_EQUAL));

    EXPECT_EQ(2, mergedBlocks.size());
    EXPECT_EQ(IndexResultBlock(0, 5), mergedBlocks.at(0));
    EXPECT_EQ(IndexResultBlock(6, 4), mergedBlocks.at(1));

    SubreadLookupData subreadIndex;
    subreadIndex.fileOffset_ = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    subreadIndex.ApplyOffsets(mergedBlocks);

    EXPECT_EQ(2,  mergedBlocks.size());
    EXPECT_EQ(0,  mergedBlocks.at(0).virtualOffset_);
    EXPECT_EQ(5,  mergedBlocks.at(0).numReads_);
    EXPECT_EQ(60, mergedBlocks.at(1).virtualOffset_);
    EXPECT_EQ(4,  mergedBlocks.at(1).numReads_);
}

TEST(PacBioIndexTest, LookupMulti)
{
    using PacBio::BAM::CompareType;
    using PacBio::BAM::IndexList;
    using PacBio::BAM::IndexResultBlock;
    using PacBio::BAM::IndexResultBlocks;
    using PacBio::BAM::SubreadField;
    using PacBio::BAM::internal::mergedIndexBlocks;
    using PacBio::BAM::internal::SubreadLookupData;
    using PacBio::BAM::internal::UnorderedLookup;

    UnorderedLookup<int32_t>::ContainerType uRawData;
    uRawData[11] = { 0, 3, 4 };
    uRawData[20] = { 1 };
    uRawData[42] = { 2, 7, 8 };
    uRawData[10] = { 5 };
    uRawData[12] = { 6 };
    uRawData[99] = { 9 };

    SubreadLookupData subreadIndex;
    subreadIndex.rgId_ = UnorderedLookup<int32_t>(uRawData);
    subreadIndex.fileOffset_ = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };

    const std::vector<int32_t> whitelist = { 11, 42, 20 };
    const auto indices = subreadIndex.IndicesMulti(SubreadField::RG_ID, whitelist);

    IndexResultBlocks mergedBlocks = mergedIndexBlocks(indices);
    subreadIndex.ApplyOffsets(mergedBlocks);

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

    // rgId == x
    const IndexResultBlocks rgResult = index.Lookup(ReadGroupIndexRequest(-1197849594));
    EXPECT_EQ(1, rgResult.size());
    EXPECT_EQ(0, rgResult.at(0).firstIndex_);
    EXPECT_EQ(4, rgResult.at(0).numReads_);
    EXPECT_EQ(35323904, rgResult.at(0).virtualOffset_);

    // rg != x
    const IndexResultBlocks notRgResult = index.Lookup(ReadGroupIndexRequest(-1197849594, CompareType::NOT_EQUAL));
    EXPECT_TRUE(notRgResult.empty());

    // tEnd <= x
    const IndexResultBlocks tEndLteResult = index.Lookup(ReferenceEndIndexRequest(9900, CompareType::LESS_THAN_EQUAL));
    EXPECT_EQ(1, tEndLteResult.size());
    EXPECT_EQ(2, tEndLteResult.at(0).firstIndex_);
    EXPECT_EQ(2, tEndLteResult.at(0).numReads_);
    EXPECT_EQ(35339466, tEndLteResult.at(0).virtualOffset_);

    // tEnd >= x
    const IndexResultBlocks tEndGteResult = index.Lookup(ReferenceEndIndexRequest(9900, CompareType::GREATER_THAN_EQUAL));
    EXPECT_EQ(2, tEndGteResult.size());
    EXPECT_EQ(0, tEndGteResult.at(0).firstIndex_);
    EXPECT_EQ(2, tEndGteResult.at(0).numReads_);
    EXPECT_EQ(35323904, tEndGteResult.at(0).virtualOffset_);
    EXPECT_EQ(3, tEndGteResult.at(1).firstIndex_);
    EXPECT_EQ(1, tEndGteResult.at(1).numReads_);
    EXPECT_EQ(35351517, tEndGteResult.at(1).virtualOffset_);

    // strand query
    const IndexResultBlocks forward = index.Lookup(StrandIndexRequest(Strand::FORWARD));
    EXPECT_EQ(2, forward.size());
    EXPECT_EQ(0, forward.at(0).firstIndex_);
    EXPECT_EQ(1, forward.at(0).numReads_);
    EXPECT_EQ(35323904, forward.at(0).virtualOffset_);
    EXPECT_EQ(2, forward.at(1).firstIndex_);
    EXPECT_EQ(1, forward.at(1).numReads_);
    EXPECT_EQ(35339466, forward.at(1).virtualOffset_);

    const IndexResultBlocks reverse = index.Lookup(StrandIndexRequest(Strand::REVERSE));
    EXPECT_EQ(2, reverse.size());
    EXPECT_EQ(1, reverse.at(0).firstIndex_);
    EXPECT_EQ(1, reverse.at(0).numReads_);
    EXPECT_EQ(35327454, reverse.at(0).virtualOffset_);
    EXPECT_EQ(3, reverse.at(1).firstIndex_);
    EXPECT_EQ(1, reverse.at(1).numReads_);
    EXPECT_EQ(35351517, reverse.at(1).virtualOffset_);

    // query data field that is not in the PBI
    const IndexResultBlocks missing = index.Lookup(BarcodeQualityIndexRequest(77, CompareType::GREATER_THAN));
    EXPECT_TRUE(missing.empty());
}

TEST(PacBioIndexTest, LookupByZmw)
{
    BamFile f(tests::Data_Dir + "/dataset/bam_mapping.bam");
    f.EnsurePacBioIndexExists();

    PbiIndex index(f.PacBioIndexFilename());

    const IndexResultBlocks blocks = index.Lookup(ZmwIndexRequest(20000, CompareType::LESS_THAN));
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
    EXPECT_EQ(32654594, block0.virtualOffset_);

    const IndexResultBlock& block1 = blocks.at(1);
    EXPECT_EQ(6, block1.firstIndex_);
    EXPECT_EQ(3, block1.numReads_);
    EXPECT_EQ(32670126, block1.virtualOffset_);

    const IndexResultBlock& block2 = blocks.at(2);
    EXPECT_EQ(10, block2.firstIndex_);
    EXPECT_EQ(3,  block2.numReads_);
    EXPECT_EQ(1389300730, block2.virtualOffset_);

    const IndexResultBlock& block3 = blocks.at(3);
    EXPECT_EQ(14, block3.firstIndex_);
    EXPECT_EQ(1,  block3.numReads_);
    EXPECT_EQ(1389323725, block3.virtualOffset_);

    const IndexResultBlock& block4 = blocks.at(4);
    EXPECT_EQ(22, block4.firstIndex_);
    EXPECT_EQ(2,  block4.numReads_);
    EXPECT_EQ(2542600192, block4.virtualOffset_);
}

TEST(PacBioIndexTest, LookupMultiZmw)
{
    BamFile f(tests::Data_Dir + "/dataset/bam_mapping.bam");
    f.EnsurePacBioIndexExists();

    PbiIndex index(f.PacBioIndexFilename());

    const std::vector<int32_t> whitelist = { 13473, 38025 };
    const ZmwIndexMultiRequest request(whitelist);
    const IndexResultBlocks& blocks = index.Lookup(request);

    EXPECT_EQ(3, blocks.size());

    const IndexResultBlock& block0 = blocks.at(0);
    EXPECT_EQ(6, block0.firstIndex_);
    EXPECT_EQ(2, block0.numReads_);
    EXPECT_EQ(32670126, block0.virtualOffset_);

    const IndexResultBlock& block1 = blocks.at(1);
    EXPECT_EQ(13, block1.firstIndex_);
    EXPECT_EQ(2, block1.numReads_);
    EXPECT_EQ(1389310464, block1.virtualOffset_);

    const IndexResultBlock& block2 = blocks.at(2);
    EXPECT_EQ(19, block2.firstIndex_);
    EXPECT_EQ(1,  block2.numReads_);
    EXPECT_EQ(1389340436, block2.virtualOffset_);
}
