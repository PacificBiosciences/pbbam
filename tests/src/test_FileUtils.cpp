// Author: Derek Barnett

#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <boost/algorithm/string.hpp>

#include "PbbamTestData.h"

#include <pbbam/../../src/FileUtils.h>
#include <pbbam/../../src/TimeUtils.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;

TEST(FileUtilsTest, ExistsOk)
{
    EXPECT_FALSE(FileUtils::Exists("does_not_exist.txt"));

    const std::string tmp = PbbamTestsConfig::GeneratedData_Dir + "/pbbam_exists_check.tmp";
    const std::string cmd = std::string("touch ") + tmp;
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
    const auto nowSeconds = std::chrono::duration_cast<std::chrono::seconds>(nowDuration).count();

    const std::string tmp = PbbamTestsConfig::GeneratedData_Dir + "/pbbam_lastmod_check.tmp";
    const std::string rmCmd = std::string("rm ") + tmp;
    const std::string touchCmd = std::string("touch  ") + tmp;
    const auto ret = system(rmCmd.c_str());
    UNUSED(ret);
    ASSERT_EQ(0, system(touchCmd.c_str()));

    const auto stamp = FileUtils::LastModified(tmp);
    const auto stampDuration = stamp.time_since_epoch();
    const auto stampSeconds =
        std::chrono::duration_cast<std::chrono::seconds>(stampDuration).count();

    EXPECT_LE(nowSeconds, stampSeconds);
}

TEST(FileUtilsTest, ResolvedFilePathOk)
{
    const std::string testFrom = "/path/to/myDir";

    // "raw" filenames - no URI scheme

    const std::string absolutePath = "/absolute/path/to/file.txt";
    const std::string relativePath = "../relative/path/to/file.txt";
    const std::string noPathFn = "file.txt";

    const std::string resolvedAbsolutePath = FileUtils::ResolvedFilePath(absolutePath, testFrom);
    const std::string resolvedRelativePath = FileUtils::ResolvedFilePath(relativePath, testFrom);
    const std::string resolvedNoPath = FileUtils::ResolvedFilePath(noPathFn, testFrom);
    const std::string resolvedAbsolutePath_defaultFrom = FileUtils::ResolvedFilePath(absolutePath);
    const std::string resolvedRelativePath_defaultFrom = FileUtils::ResolvedFilePath(relativePath);
    const std::string resolvedNoPath_defaultFrom = FileUtils::ResolvedFilePath(noPathFn);

    EXPECT_EQ("/absolute/path/to/file.txt", resolvedAbsolutePath);
    EXPECT_EQ("/path/to/myDir/../relative/path/to/file.txt", resolvedRelativePath);
    EXPECT_EQ("/path/to/myDir/file.txt", resolvedNoPath);

    EXPECT_EQ("/absolute/path/to/file.txt", resolvedAbsolutePath_defaultFrom);
    EXPECT_EQ("./../relative/path/to/file.txt", resolvedRelativePath_defaultFrom);
    EXPECT_EQ("./file.txt", resolvedNoPath_defaultFrom);

    // filenames with URI scheme ("file://")

    const std::string absoluteSchemeFn = "file:///absolute/path/to/file.txt";
    const std::string relativeSchemeFn = "file://../relative/path/to/file.txt";
    const std::string noPathSchemeFn = "file://file.txt";

    const std::string resolvedAbsoluteSchemePath =
        FileUtils::ResolvedFilePath(absoluteSchemeFn, testFrom);
    const std::string resolvedRelativeSchemePath =
        FileUtils::ResolvedFilePath(relativeSchemeFn, testFrom);
    const std::string resolvedNoPathSchemeFn =
        FileUtils::ResolvedFilePath(noPathSchemeFn, testFrom);
    const std::string resolvedAbsoluteSchemePath_defaultFrom =
        FileUtils::ResolvedFilePath(absoluteSchemeFn);
    const std::string resolvedRelativeSchemePath_defaultFrom =
        FileUtils::ResolvedFilePath(relativeSchemeFn);
    const std::string resolvedNoPathSchemeFn_defaultFrom =
        FileUtils::ResolvedFilePath(noPathSchemeFn);

    EXPECT_EQ("/absolute/path/to/file.txt", resolvedAbsoluteSchemePath);
    EXPECT_EQ("/path/to/myDir/../relative/path/to/file.txt", resolvedRelativeSchemePath);
    EXPECT_EQ("/path/to/myDir/file.txt", resolvedNoPathSchemeFn);

    EXPECT_EQ("/absolute/path/to/file.txt", resolvedAbsoluteSchemePath_defaultFrom);
    EXPECT_EQ("./../relative/path/to/file.txt", resolvedRelativeSchemePath_defaultFrom);
    EXPECT_EQ("./file.txt", resolvedNoPathSchemeFn_defaultFrom);
}

TEST(FileUtilsTest, SizeOk)
{
    const std::string tmp = PbbamTestsConfig::GeneratedData_Dir + "/pbbam_empty_file.tmp";
    const std::string cmd = std::string("touch ") + tmp;
    ASSERT_EQ(0, system(cmd.c_str()));
    EXPECT_EQ(0, FileUtils::Size(tmp));

    EXPECT_THROW(FileUtils::Size("does_not_exist.txt"), std::runtime_error);
}

// ####################################################################################################
// The code below is part of a simple check whether or not a (Windows-only) file path is absolute.
//
// NOTE: (and this is admittedly brittle for maintenance, but) the internal methods used are literally
// copied here for direct driving. There's likely a better way going forward, than the manual copy/paste.
// But in the absence of a similar runtime environment to build in & test against, while
// the motivating behavior is blocking other work, this lets me get the fix in their hands ASAP and still
// have some test code poking it beforehand. -DB
//
namespace test_windows {

static std::string removeFileUriScheme(const std::string& uri)
{
    assert(!uri.empty());

    auto schemeLess = uri;
    const auto fileScheme = std::string{"file://"};
    const auto schemeFound = schemeLess.find(fileScheme);
    if (schemeFound != std::string::npos) {
        if (schemeFound != 0) throw std::runtime_error("Malformed URI: scheme not at beginning");
        schemeLess = schemeLess.substr(fileScheme.size());
    }
    return schemeLess;
}

static std::string removeDiskName(const std::string& filePath)
{
    if (filePath.size() >= 2) {
        const char firstChar = filePath.at(0);
        if ((isalpha(firstChar) != 0) && (filePath.at(1) == ':')) return filePath.substr(2);
    }
    return filePath;
}

static const char native_pathSeparator = '\\';

static bool native_pathIsAbsolute(const std::string& filePath)
{
    assert(!filePath.empty());

    // if starts with single slash or double slash [cases 1,3]
    if (boost::algorithm::starts_with(filePath, "\\")) return true;

    // if starts with single or double-dots -> not absolute [case 4 + ".\file.txt"]
    if (boost::algorithm::starts_with(filePath, ".")) return false;

    // if starts with drive name and colon ("C:\foo\bar.txt")
    if (filePath.size() >= 2) {
        const char firstChar = filePath.at(0);
        if ((isalpha(firstChar) != 0) && (filePath.at(1) == ':'))
            return native_pathIsAbsolute(removeDiskName(filePath));
    }

    // otherwise, likely relative
    return false;
}

static std::string native_resolvedFilePath(const std::string& filePath, const std::string& from)
{
    // strip file:// scheme if present
    auto schemeLess = removeFileUriScheme(filePath);

    // if empty or already absolute path, just return it
    // upfront empty check simplifies further parsing logic
    if (schemeLess.empty() || native_pathIsAbsolute(schemeLess)) return schemeLess;

    // else make relative from the provided 'from' directory
    //
    // first pop disk name, then any leading single-dot '.'
    //
    // since we're prepending the 'from' directory, we can remove
    // any leading './' form our file path. this may just mean that
    // we pop it off to add it right back (when from == '.'), but this
    // keeps it consistent with other 'from' parent directories
    //
    schemeLess = removeDiskName(schemeLess);

    const bool thisDirAtStart = (schemeLess.find(".") == 0);
    if (thisDirAtStart) {
        if (schemeLess.find(native_pathSeparator) == 1) schemeLess = schemeLess.substr(2);
    }
    return from + native_pathSeparator + schemeLess;
}

}  // namespace test_windows

TEST(FileUtilsTest, WindowsPathsOk)
{
    {  // remove disk name

        // "C:\tmp.txt"
        std::string f1 = "C:\\tmp.txt";
        EXPECT_EQ(std::string("\\tmp.txt"), test_windows::removeDiskName(f1));

        // "C:tmp.txt"
        std::string f2 = "C:tmp.txt";
        EXPECT_EQ(std::string("tmp.txt"), test_windows::removeDiskName(f2));

        // "\tmp.txt"
        std::string f3 = "\\tmp.txt";
        EXPECT_EQ(f3, test_windows::removeDiskName(f3));

        // "tmp.txt"
        std::string f4 = "tmp.txt";
        EXPECT_EQ(f4, test_windows::removeDiskName(f4));
    }

    {  // isAbsolute ?

        // "\\server\path\to\tmp.txt"
        EXPECT_TRUE(test_windows::native_pathIsAbsolute("\\\\server\\path\\to\tmp.txt"));

        // "..\tmp.txt"
        EXPECT_FALSE(test_windows::native_pathIsAbsolute("..\\tmp.txt"));

        // ".\tmp.txt"
        EXPECT_FALSE(test_windows::native_pathIsAbsolute(".\\tmp.txt"));

        // "C:\path\to\tmp.txt"
        EXPECT_TRUE(test_windows::native_pathIsAbsolute("C:\\path\\to\\tmp.txt"));

        // "C:..\path\to\tmp.txt"
        EXPECT_FALSE(test_windows::native_pathIsAbsolute("C:..\\path\\to\\tmp.txt"));
    }

    {  // resolve file path

        const std::string myRootDir = "C:\\path\\to\\myRootDir";

        // "\\server\path\to\tmp.txt"
        const std::string fn1 = "\\\\server\\path\\to\tmp.txt";
        const std::string fn1_expected = fn1;
        EXPECT_EQ(fn1_expected, test_windows::native_resolvedFilePath(fn1, myRootDir));

        // "..\tmp.txt"
        const std::string fn2 = "..\\tmp.txt";
        const std::string fn2_expected = "C:\\path\\to\\myRootDir\\..\\tmp.txt";
        EXPECT_EQ(fn2_expected, test_windows::native_resolvedFilePath(fn2, myRootDir));

        // ".\tmp.txt"
        const std::string fn3 = ".\\tmp.txt";
        const std::string fn3_expected = "C:\\path\\to\\myRootDir\\tmp.txt";
        EXPECT_EQ(fn3_expected, test_windows::native_resolvedFilePath(fn3, myRootDir));

        // "C:\path\to\tmp.txt"
        const std::string fn4 = "C:\\path\\to\\tmp.txt";
        const std::string fn4_expected = fn4;
        EXPECT_EQ(fn4_expected, test_windows::native_resolvedFilePath(fn4, myRootDir));

        // "C:..\path\to\tmp.txt"
        const std::string fn5 = "C:..\\path\\to\\tmp.txt";
        const std::string fn5_expected = "C:\\path\\to\\myRootDir\\..\\path\\to\\tmp.txt";
        EXPECT_EQ(fn5_expected, test_windows::native_resolvedFilePath(fn5, myRootDir));

        // "C:tmp.txt"
        const std::string fn6 = "C:tmp.txt";
        const std::string fn6_expected = "C:\\path\\to\\myRootDir\\tmp.txt";
        EXPECT_EQ(fn6_expected, test_windows::native_resolvedFilePath(fn6, myRootDir));
        EXPECT_EQ(fn3_expected,
                  test_windows::native_resolvedFilePath(
                      fn6, myRootDir));  // our path is equivalent to fn3's "./temp.txt"
    }
}
//
// ####################################################################################################
