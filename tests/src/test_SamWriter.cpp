// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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

#include <gtest/gtest.h>
#include <pbbam/SamWriter.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include "PbbamTestData.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(SamWriterTest, HeaderOk)
{
    // setup header
    const string hdrText = {
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.3\n"
        "@RG\tID:6002b307\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
        "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
        "PU:test\tPM:SEQUEL\n"};

    EXPECT_NO_THROW({
        // write header to file
        const string generatedFn = PbbamTestsConfig::GeneratedData_Dir + "/samwriter_hdr_only.sam";
        {
            const BamHeader inputHeader(hdrText);
            SamWriter writer(generatedFn, inputHeader);
            //            ()writer;
        };

        // check header
        {
            ifstream f(generatedFn);
            const string text((istreambuf_iterator<char>(f)), istreambuf_iterator<char>());
            EXPECT_EQ(hdrText, text);
        }

        // clean up
        remove(generatedFn.c_str());
    });
}

TEST(SamWriterTest, SingleRecordOk)
{

    // setup header
    const string hdrLine1 = {"@HD\tVN:1.1\tSO:unknown\tpb:3.0.3"};
    const string hdrLine2 = {
        "@RG\tID:6002b307\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
        "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
        "PU:test\tPM:SEQUEL"};
    const string hdrText = hdrLine1 + "\n" + hdrLine2 + "\n";
    const BamHeader inputHeader(hdrText);

    // setup record
    BamRecord record(inputHeader);
    record.Impl().Name("test/100/0_5");
    record.Impl().SetSequenceAndQualities("ACGTC", 5, "@@@@@");
    record.Impl().CigarData("");
    record.Impl().Bin(0);
    record.Impl().Flag(0);
    record.Impl().InsertSize(0);
    record.Impl().MapQuality(0);
    record.Impl().MatePosition(-1);
    record.Impl().MateReferenceId(-1);
    record.Impl().Position(-1);
    record.Impl().ReferenceId(-1);
    record.Impl().SetMapped(false);

    TagCollection tags;
    tags["zm"] = static_cast<int32_t>(100);
    tags["qs"] = static_cast<Position>(0);
    tags["qe"] = static_cast<Position>(5);
    tags["np"] = static_cast<int32_t>(1);
    tags["rq"] = static_cast<float>(0.6);
    tags["RG"] = std::string{"6002b307"};
    tags["sn"] = vector<float>{0.2f, 0.2f, 0.2f, 0.2f};
    record.Impl().Tags(tags);

    const string expectedSamRecord = {
        "test/100/0_5\t4\t*\t0\t0\t*\t*\t0\t0\tACGTC\t@@@@@\tRG:Z:6002b307\t"
        "np:i:1\tqe:i:5\tqs:i:0\trq:f:0.6\tsn:B:f,0.2,0.2,0.2,0.2\tzm:i:100"};

    EXPECT_NO_THROW({
        // write data to file
        const string generatedFn =
            PbbamTestsConfig::GeneratedData_Dir + "/samwriter_hdr_and_record.sam";
        {
            SamWriter writer(generatedFn, inputHeader);
            writer.Write(record);
        };

        // check header & record
        {
            ifstream f(generatedFn);
            string line1;
            string line2;
            string line3;
            std::getline(f, line1);
            std::getline(f, line2);
            std::getline(f, line3);
            EXPECT_EQ(hdrLine1, line1);
            EXPECT_EQ(hdrLine2, line2);
            EXPECT_EQ(expectedSamRecord, line3);
        }

        // cleanup
        remove(generatedFn.c_str());
    });
}
