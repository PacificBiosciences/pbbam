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
#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(BamWriterTest, SingleWrite_UserRecord)
{
    const string fullName = "test/100/0_5";
    const string rgId     = "6002b307";
    const vector<float> expectedSnr = {0.2,0.2,0.2,0.2};

    // setup header
    const string hdrText = {
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:6002b307\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
             "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
             "PU:test\tPM:SEQUEL\n"
    };
    BamHeader inputHeader(hdrText);

    // setup record
    BamRecord bamRecord(inputHeader);
    bamRecord.impl_.Name(fullName);
    bamRecord.impl_.SetSequenceAndQualities("ACGTC", 5);
    bamRecord.impl_.CigarData("");
    bamRecord.impl_.Bin(0);
    bamRecord.impl_.Flag(0);
    bamRecord.impl_.InsertSize(0);
    bamRecord.impl_.MapQuality(0);
    bamRecord.impl_.MatePosition(-1);
    bamRecord.impl_.MateReferenceId(-1);
    bamRecord.impl_.Position(-1);
    bamRecord.impl_.ReferenceId(-1);
    bamRecord.impl_.SetMapped(false);

    TagCollection tags;
    tags["zm"] = static_cast<int32_t>(100);
    tags["qs"] = static_cast<Position>(0);
    tags["qe"] = static_cast<Position>(5);
    tags["np"] = static_cast<int32_t>(1);
    tags["rq"] = static_cast<float>(0.6);
    tags["RG"] = rgId;
    tags["sn"] = expectedSnr;
    bamRecord.impl_.Tags(tags);

    // write record to file
    const string generatedBamFn = tests::GeneratedData_Dir + "/bamwriter_generated.bam";
    {
        BamWriter writer(generatedBamFn, inputHeader);
        writer.Write(bamRecord);
    }

    // check written header
    BamFile file(generatedBamFn);
    const auto header = file.Header();
    EXPECT_EQ(std::string("1.1"),     header.Version());
    EXPECT_EQ(std::string("unknown"), header.SortOrder());
    EXPECT_EQ(std::string("3.0.1"),   header.PacBioBamVersion());

    // check written record
    EntireFileQuery entireFile(file);
    auto firstIter = entireFile.begin();
    auto record = *firstIter;
    EXPECT_EQ(std::string("ACGTC"),        record.Sequence());
    EXPECT_EQ(std::string("test/100/0_5"), record.FullName());
    EXPECT_TRUE(record.HasHoleNumber());
    EXPECT_TRUE(record.HasNumPasses());
    EXPECT_TRUE(record.HasQueryEnd());
    EXPECT_TRUE(record.HasQueryStart());
    EXPECT_TRUE(record.HasReadAccuracy());
    EXPECT_TRUE(record.HasSignalToNoise());
    EXPECT_EQ(100, record.HoleNumber());
    EXPECT_EQ(1,   record.NumPasses());
    EXPECT_EQ(0,   record.QueryStart());
    EXPECT_EQ(5,   record.QueryEnd());
    EXPECT_EQ(expectedSnr, record.SignalToNoise());
    EXPECT_EQ(rgId, record.ReadGroupId());

    // clean up
    remove(generatedBamFn.c_str());
}
