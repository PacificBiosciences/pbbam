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
    // active development - just throws for now
    BamFile bamFile(test2BamFn);
    EXPECT_THROW(PbiIndex index(bamFile.PacBioIndexFilename()), std::exception);


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
