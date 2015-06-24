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
#include <htslib/sam.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <chrono>
#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

// put any BamWriter-only API tests here (error handling, etc.)
//
// plain ol' read & dump is in test_EndToEnd.cpp

const string generatedBamFn = tests::Data_Dir + "/generated.bam";

struct ResultPacket
{
    std::string name;
    char*       bases;
    char*       overallQv;
    size_t      length;
    int         zmwNum;
    int         startPos;
    BamRecord   bamRecord;

    ResultPacket() = default;

    ResultPacket(ResultPacket&& src)
    {
        name            = std::move(src.name);

        bases           = src.bases;
        overallQv       = src.overallQv;
        length          = src.length;

        zmwNum          = src.zmwNum;
        startPos        = src.startPos;

        src.bases     = 0;
        src.overallQv = 0;

        bamRecord = std::move(src.bamRecord);
    }
    // Copy constructor
    ResultPacket(const ResultPacket&) = delete;
    // Move assignment constructor
    ResultPacket& operator=(ResultPacket&&) = delete;
    // Copy assignment constructor
    ResultPacket& operator=(const ResultPacket&) = delete;
    // Destructor
    ~ResultPacket()
    {
        // delete [] basesBam;
        if (bases != 0) delete [] bases;
        if (overallQv != 0) delete [] overallQv;
    }
};

TEST(BamWriterTest, SingleWrite_UserRecord)
{
    //Writing a ResultPacket in Workflow.h:
    ResultPacket result;
    result.zmwNum = 42;
    result.name = "ZMW\\"+std::to_string(42);
    auto length = 5;

    result.bases     = (char*) calloc(length,1);
    result.overallQv = (char*) calloc(length,1);
    // FILL WITH CONTENT
    result.bases[0] = 'A';
    result.bases[1] = 'C';
    result.bases[2] = 'G';
    result.bases[3] = 'T';
    result.bases[4] = 'C';
    result.overallQv[0] = ']';
    result.overallQv[1] = '6';
    result.overallQv[2] = '4';
    result.overallQv[3] = '@';
    result.overallQv[4] = '<';

    // Encode data to BamAlignment
    result.bamRecord.impl_.Name(result.name);
    result.bamRecord.impl_.SetSequenceAndQualities(result.bases, length);
    result.bamRecord.impl_.CigarData("");
    result.bamRecord.impl_.Bin(0);
    result.bamRecord.impl_.Flag(0);
    result.bamRecord.impl_.InsertSize(0);
    result.bamRecord.impl_.MapQuality(0);
    result.bamRecord.impl_.MatePosition(-1);
    result.bamRecord.impl_.MateReferenceId(-1);
    result.bamRecord.impl_.Position(-1);
    result.bamRecord.impl_.ReferenceId(-1);

    std::vector<uint8_t> subQv = std::vector<uint8_t>({34, 5, 125});

    TagCollection tags;
    tags["SQ"] = subQv;

    Tag asciiTag('J');
    asciiTag.Modifier(TagModifier::ASCII_CHAR);

    // add ASCII tag via TagCollection
    tags["a1"] = asciiTag;
    result.bamRecord.impl_.Tags(tags);

    // add ASCII tag via BamRecordImpl
    Tag asciiTag2('K');
    asciiTag2.Modifier(TagModifier::ASCII_CHAR);
    result.bamRecord.impl_.AddTag("a2", asciiTag2);

    BamHeader headerSubreads;
    headerSubreads.Version("1.1")
                  .SortOrder("coordinate");

    EXPECT_NO_THROW ({
        BamWriter writer(generatedBamFn, headerSubreads);
        writer.Write(result.bamRecord);
    });

    EXPECT_NO_THROW ({
        BamFile file(generatedBamFn);
        EXPECT_EQ(std::string("1.1"),        file.Header().Version());
        EXPECT_EQ(std::string("coordinate"), file.Header().SortOrder());

        EntireFileQuery entireFile(file);
        for (const BamRecord& record : entireFile) {
            const BamRecordImpl& impl = record.Impl();

            EXPECT_EQ(std::string("ACGTC"),   impl.Sequence());
            EXPECT_EQ(std::string("ZMW\\42"), impl.Name());

            const TagCollection& implTags = impl.Tags();
            EXPECT_TRUE(implTags.Contains("SQ"));
            EXPECT_TRUE(implTags.Contains("a1"));
            EXPECT_TRUE(implTags.Contains("a2"));

            const Tag sqTag = impl.TagValue("SQ");
            const Tag a1Tag = impl.TagValue("a1");
            const Tag a2Tag = impl.TagValue("a2");
            EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), sqTag.ToUInt8Array());
            EXPECT_EQ('J', a1Tag.ToAscii());
            EXPECT_EQ('K', a2Tag.ToAscii());

            // just check first record
            break;
        }
    });

    remove(generatedBamFn.c_str());
}

//#define SEQ_LENGTH  7000
//#define NUM_RECORDS 100000

//static const std::string& TEST_SEQUENCE  = std::string(SEQ_LENGTH, 'G');
//static const std::string& TEST_QUALITIES = std::string(SEQ_LENGTH, '=');
//static const std::string& TEST_NAME      = std::string(SEQ_LENGTH, '/');
//static const std::string& TEST_TAGDATA   = std::string(SEQ_LENGTH, '2');

//TEST(BamWriterTest, CheckTiming)
//{
//    TagCollection tags;
//    tags["aa"] = TEST_TAGDATA;
//    tags["bb"] = TEST_TAGDATA;
//    tags["cc"] = TEST_TAGDATA;
//    tags["dd"] = TEST_TAGDATA;
//    tags["ee"] = TEST_TAGDATA;
//    tags["ff"] = TEST_TAGDATA;

//    BamRecord record;
//    record.SetSequenceAndQualities(TEST_SEQUENCE, TEST_QUALITIES);
//    record.Name(TEST_NAME);
//    record.Tags(tags);

//    SamHeader header;
//    header.pacbioBamVersion = "3.0b3";

//    BamWriter writer("fake.bam", header);
//    if(!writer) {
//        std::cout << "BamWriter::WriteError: " << writer.Error() << endl;
//    }


//    auto start = std::chrono::steady_clock::now();
//    for (size_t i = 0; i < NUM_RECORDS; ++i)
//        writer.Write(record);
//    writer.Close();
//    auto end = std::chrono::steady_clock::now();

//    auto diff = end - start;
//    std::cout << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
//}

