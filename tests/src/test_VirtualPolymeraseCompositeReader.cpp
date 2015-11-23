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
#include <pbbam/virtual/VirtualPolymeraseCompositeReader.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace tests {

static
size_t NumVirtualRecords(const string& primaryBamFn,
                         const string& scrapsBamFn)
{
    VirtualPolymeraseReader reader(primaryBamFn, scrapsBamFn);
    size_t count = 0;
    while (reader.HasNext()) {
        const auto record = reader.Next();
        (void)record;
        ++count;
    }
    return count;
}

} // namespace tests
} // namespace BAM
} // namespace PacBio

TEST(VirtualPolymeraseCompositeReaderTest, DataSetOk)
{
    // dataset contains these resources (subreads/scraps + hqregion/scraps BAMs)
    const string primaryFn1 = tests::Data_Dir + "/polymerase/production.subreads.bam";
    const string scrapsFn1  = tests::Data_Dir + "/polymerase/production.scraps.bam";
    const string primaryFn2 = tests::Data_Dir + "/polymerase/production_hq.hqregion.bam";
    const string scrapsFn2  = tests::Data_Dir + "/polymerase/production_hq.scraps.bam";
    const size_t numExpectedRecords =
            tests::NumVirtualRecords(primaryFn1, scrapsFn1) +
            tests::NumVirtualRecords(primaryFn2, scrapsFn2);

    const string datasetFn = tests::Data_Dir +
            "/polymerase/multiple_resources.subread.dataset.xml";

    DataSet ds(datasetFn);
    VirtualPolymeraseCompositeReader reader(ds);
    size_t numObservedRecords = 0;
    while (reader.HasNext()) {
        const auto record = reader.Next();
        (void)record;
        ++numObservedRecords;
    }
    EXPECT_EQ(numExpectedRecords, numObservedRecords);
}

TEST(VirtualPolymeraseCompositeReaderTest, EmptyDataSetOk)
{
    VirtualPolymeraseCompositeReader reader(DataSet{});
    EXPECT_FALSE(reader.HasNext());
}
