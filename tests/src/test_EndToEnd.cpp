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

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <gtest/gtest.h>
#include <htslib/sam.h>

#define private public
#define protected public

#include "PbbamTestData.h"

#include <pbbam/BamFile.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace EndToEndTests {

struct Bam1Deleter
{
    void operator()(bam1_t* b) {
        if (b)
            bam_destroy1(b);
        b = nullptr;
    }
};

struct SamFileDeleter
{
    void operator()(samFile* file) {
        if (file)
            sam_close(file);
        file = nullptr;
    }
};

struct BamHdrDeleter
{
    void operator()(bam_hdr_t* hdr) {
        if (hdr)
            bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

const string inputBamFn        = PbbamTestsConfig::Data_Dir + "/aligned.bam";
const string goldStandardSamFn = PbbamTestsConfig::Data_Dir + "/aligned.sam";
const string generatedBamFn    = PbbamTestsConfig::GeneratedData_Dir + "/generated.bam";
const string generatedSamFn    = PbbamTestsConfig::GeneratedData_Dir + "/generated.sam";
const vector<string> generatedFiles = { generatedBamFn, generatedSamFn };

static inline
int RunBam2Sam(const string& bamFn,
               const string& samFn,
               const string& args = string())
{
    stringstream s;
    s << PbbamTestsConfig::Bam2Sam << " " << args << " " << bamFn << " > " << samFn;
    return system(s.str().c_str());
}

static inline
int RunDiff(const string& fn1, const string& fn2)
{
    stringstream s;
    s << "diff " << fn1 << " " << fn2;
    return system(s.str().c_str());
}

static inline
void Remove(const vector<string>& files)
{
    for (const auto& fn : files)
        remove(fn.c_str());
}

static inline
void CheckGeneratedOutput(void)
{
    // convert to sam & diff against gold standard
    const int convertRet = RunBam2Sam(generatedBamFn, generatedSamFn);
    const int diffRet    = RunDiff(goldStandardSamFn, generatedSamFn);
    EXPECT_EQ(0, convertRet);
    EXPECT_EQ(0, diffRet);

    // clean up
    Remove(generatedFiles);
}

} // namespace EndToEndTests

// sanity check for rest of tests below
TEST(EndToEndTest, ReadAndWrite_PureHtslib)
{
    { // scoped to force flush & close before conversion/diff

        // open files

        unique_ptr<samFile, EndToEndTests::SamFileDeleter> inWrapper(sam_open(EndToEndTests::inputBamFn.c_str(), "r"));
        samFile* in = inWrapper.get();
        ASSERT_TRUE(in);

        unique_ptr<samFile, EndToEndTests::SamFileDeleter> outWrapper(sam_open(EndToEndTests::generatedBamFn.c_str(), "wb"));
        samFile* out = outWrapper.get();
        ASSERT_TRUE(out);

        // fetch & write header

        unique_ptr<bam_hdr_t, EndToEndTests::BamHdrDeleter> headerWrapper(sam_hdr_read(in));
        bam_hdr_t* hdr = headerWrapper.get();
        ASSERT_TRUE(hdr);
        ASSERT_EQ(0, sam_hdr_write(out, hdr));

        // fetch & write records

        unique_ptr<bam1_t, EndToEndTests::Bam1Deleter> record(bam_init1());
        bam1_t* b = record.get();
        ASSERT_TRUE(b);

        while (sam_read1(in, hdr, b) >= 0)
            sam_write1(out, hdr, b);
    }

    EndToEndTests::CheckGeneratedOutput();
}

TEST(EndToEndTest, ReadAndWrite_SingleThread)
{
    EXPECT_NO_THROW(
    {
        // open input BAM file
        BamFile bamFile(EndToEndTests::inputBamFn);

        // open output BAM file
        BamWriter writer(EndToEndTests::generatedBamFn, bamFile.Header(), BamWriter::DefaultCompression, 1);

        // copy BAM file
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile)
            writer.Write(record);
    });

    EndToEndTests::CheckGeneratedOutput();
}

TEST(EndToEndTest, ReadAndWrite_APIDefaultThreadCount)
{
    EXPECT_NO_THROW(
    {
        // open input BAM file
        BamFile bamFile(EndToEndTests::inputBamFn);

        // open output BAM file
        BamWriter writer(EndToEndTests::generatedBamFn, bamFile.Header());

        // copy BAM file
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile)
            writer.Write(record);
    });

    EndToEndTests::CheckGeneratedOutput();
}

TEST(EndToEndTest, ReadAndWrite_SystemDefaultThreadCount)
{
    EXPECT_NO_THROW(
    {
        // open input BAM file
        BamFile bamFile(EndToEndTests::inputBamFn);

        // open output BAM file
        BamWriter writer(EndToEndTests::generatedBamFn,
                         bamFile.Header(),
                         BamWriter::DefaultCompression,
                         0);

        // copy BAM file
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile)
            writer.Write(record);
    });

    EndToEndTests::CheckGeneratedOutput();
}

TEST(EndToEndTest, ReadAndWrite_UserThreadCount)
{
    EXPECT_NO_THROW(
    {
        // open input BAM file
        BamFile bamFile(EndToEndTests::inputBamFn);

        // open output BAM file
        BamWriter writer(EndToEndTests::generatedBamFn,
                         bamFile.Header(),
                         BamWriter::DefaultCompression,
                         3);

        // copy BAM file
        EntireFileQuery entireFile(bamFile);
        for (const BamRecord& record : entireFile)
            writer.Write(record);
    });

    EndToEndTests::CheckGeneratedOutput();
}
