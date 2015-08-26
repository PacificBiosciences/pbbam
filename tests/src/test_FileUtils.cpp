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

#include <gtest/gtest.h>
#include <pbbam/../../src/FileUtils.h>
#include <pbbam/../../src/TimeUtils.h>
#include <chrono>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include <iostream>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

TEST(FileUtilsTest, ExistsOk)
{
    EXPECT_FALSE(FileUtils::Exists("does_not_exist.txt"));

    const string tmp = "/tmp/pbbam_exists_check.tmp";
    const string cmd = string("touch ") + tmp;
    ASSERT_EQ(0, system(cmd.c_str()));
    EXPECT_TRUE(FileUtils::Exists(tmp));
}

TEST(FileUtilsTest, LastModifiedOk)
{
    // a little tricky to check without going a full 'mock' filesystem route, but we can approximate
    //
    // also, I can't seem to get better than second resolution (on OSX 10.9/clang at least, st_mtimespec.tv_nsec is always zero)

    const auto now = CurrentTime();
    const auto nowDuration = now.time_since_epoch();
    const auto nowSeconds = chrono::duration_cast<chrono::seconds>(nowDuration).count();

    const string tmp = "/tmp/pbbam_lastmod_check.tmp";
    const string rmCmd = string("rm ") + tmp;
    const string touchCmd = string("touch  ") + tmp;
    int ret =  system(rmCmd.c_str());
    (void)ret; // unused
    ASSERT_EQ(0, system(touchCmd.c_str()));

    const auto stamp = FileUtils::LastModified(tmp);
    const auto stampDuration = stamp.time_since_epoch();
    const auto stampSeconds = chrono::duration_cast<chrono::seconds>(stampDuration).count();

    EXPECT_LE(nowSeconds, stampSeconds);
}

TEST(FileUtilsTest, ResolvedFilePathOk)
{
    const string absolutePath = "/absolute/path/to/file.txt";
    const string relativePath = "../relative/path/to/file.txt";
    const string noPathFn     = "file.txt";

    const string testFrom = "/path/to/myDir";

    const string resolvedAbsolutePath = FileUtils::ResolvedFilePath(absolutePath, testFrom);
    const string resolvedRelativePath = FileUtils::ResolvedFilePath(relativePath, testFrom);
    const string resolvedNoPath       = FileUtils::ResolvedFilePath(noPathFn, testFrom);
    const string resolvedAbsolutePath_defaultFrom = FileUtils::ResolvedFilePath(absolutePath);
    const string resolvedRelativePath_defaultFrom = FileUtils::ResolvedFilePath(relativePath);
    const string resolvedNoPath_defaultFrom       = FileUtils::ResolvedFilePath(noPathFn);

    EXPECT_EQ("/absolute/path/to/file.txt",                  resolvedAbsolutePath);
    EXPECT_EQ("/path/to/myDir/../relative/path/to/file.txt", resolvedRelativePath);
    EXPECT_EQ("/path/to/myDir/file.txt",                     resolvedNoPath);

    EXPECT_EQ("/absolute/path/to/file.txt",     resolvedAbsolutePath_defaultFrom);
    EXPECT_EQ("./../relative/path/to/file.txt", resolvedRelativePath_defaultFrom);
    EXPECT_EQ("./file.txt",                     resolvedNoPath_defaultFrom);
}

TEST(FileUtilsTest, SizeOk)
{
    const string tmp = "/tmp/pbbam_empty_file.tmp";
    const string cmd = string("touch ") + tmp;
    ASSERT_EQ(0, system(cmd.c_str()));
    EXPECT_EQ(0, FileUtils::Size(tmp));

    EXPECT_THROW(FileUtils::Size("does_not_exist.txt"), std::runtime_error);
}
