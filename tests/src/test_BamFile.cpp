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
#include <pbbam/EntireFileQuery.h>
#include <pbbam/../../src/FileUtils.h>
#include <stdexcept>
#include <cstdlib>
#include <unistd.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(BamFileTest, NonExistentFileThrows)
{
    EXPECT_THROW(
    {
       BamFile file("does_not_exist.bam");
       (void)file;
    },
    std::exception);
}

TEST(BamFileTest, NonBamFileThrows)
{
    EXPECT_THROW(
    {
        const std::string& fn = tests::Data_Dir + "/lambdaNEB.fa.fai";
        BamFile file(fn);
        (void)file;
    },
    std::exception);
}

TEST(BamFileTest, RelativePathBamOk)
{
    const string cwd = internal::FileUtils::CurrentWorkingDirectory();
    ASSERT_EQ(0, chdir(tests::Data_Dir.c_str()));
    ASSERT_EQ(0, chdir("relative/a"));

    { // direct BAM
        BamFile file("../b/test1.bam");
        EntireFileQuery entireFile(file);
        int count = 0;
        for (const BamRecord& r : entireFile) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(10, count);
    }

    { // dataset from BAM filename
        DataSet ds("../b/test1.bam");
        EntireFileQuery entireFile(ds);
        int count = 0;
        for (const BamRecord& r : entireFile) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(10, count);
    }

    { // dataset from BamFile object
        BamFile file("../b/test1.bam");
        DataSet ds(file);
        EntireFileQuery entireFile(ds);
        int count = 0;
        for (const BamRecord& r : entireFile) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(10, count);
    }

    ASSERT_EQ(0, chdir(cwd.c_str()));
}

TEST(BamFileTest, RelativePathXmlOk)
{
    const string cwd = internal::FileUtils::CurrentWorkingDirectory();

    ASSERT_EQ(0, chdir(tests::Data_Dir.c_str()));

    {
        DataSet ds("relative/relative.xml");
        EntireFileQuery entireFile(ds);
        int count = 0;
        for (const BamRecord& r : entireFile) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(30, count);
    }

    ASSERT_EQ(0, chdir(cwd.c_str()));
}

TEST(BamFileTest, RelativePathFofnOk)
{
    const string cwd = internal::FileUtils::CurrentWorkingDirectory();
    ASSERT_EQ(0, chdir(tests::Data_Dir.c_str()));

    { // FOFN containing BAMs in different subdirs

        DataSet ds("relative/relative.fofn");
        EntireFileQuery entireFile(ds);
        int count = 0;
        for (const BamRecord& r : entireFile) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(30, count);
    }

    // NOTE: doesn't yet support a FOFN containing an XML with relative paths

//    { // FOFN containing subdir BAMs + relative.xml

//        DataSet ds("relative/relative2.fofn");
//        EntireFileQuery entireFile(ds);
//        int count = 0;
//        for (const BamRecord& r : entireFile) {
//            (void)r;
//            ++count;
//        }
//        EXPECT_EQ(60, count);
//    }

    ASSERT_EQ(0, chdir(cwd.c_str()));
}
