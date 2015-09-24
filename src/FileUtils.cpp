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

#include "FileUtils.h"
#include "StringUtils.h"
#include <boost/algorithm/string.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <cassert>
#include <sys/stat.h>
#include <unistd.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// pops "file://" scheme off the front of a URI/filepath, if found
static string removeFileUriScheme(const string& uri)
{
    assert(!uri.empty());

    auto schemeLess = uri;
    const auto fileScheme = string{"file://"};
    const auto schemeFound = schemeLess.find(fileScheme);
    if (schemeFound != string::npos) {
        if (schemeFound != 0)
            throw runtime_error("Malformed URI: scheme not at beginning");
        schemeLess = schemeLess.substr(fileScheme.size());
    }
    return schemeLess;
}

#ifdef PBBAM_WIN_FILEPATHS

static
string removeDiskName(const string& filePath)
{
    if (filePath.size() >= 2) {
        const char firstChar = filePath.at(0);
        if ((isalpha(firstChar) != 0) && (filePath.at(1) == ':'))
            return filePath.substr(2);
    }
    return filePath;
}

static const char native_pathSeparator = '\\';

static bool native_pathIsAbsolute(const string& filePath)
{
    assert(!filePath.empty());

    // if starts with single slash or double slash
    if (boost::algorithm::starts_with(filePath, "\\"))
        return true;

    // if starts with single or double-dots -> not absolute
    if (boost::algorithm::starts_with(filePath, "."))
        return false;

    // if starts with disk drive name and colon ("C:\foo\bar.txt")
    // strip the drive name and check to see if the remaining path is absolute
    if (filePath.size() >= 2) {
        const char firstChar = filePath.at(0);
        if ((isalpha(firstChar) != 0) && (filePath.at(1) == ':'))
            return native_pathIsAbsolute(removeDiskName(filePath));
    }

    // otherwise, likely relative
    return false;
}

static string native_resolvedFilePath(const string& filePath,
                                      const string& from)
{
    // strip file:// scheme if present
    auto schemeLess = removeFileUriScheme(filePath);

    // if empty or already absolute path, just return it
    // upfront empty check simplifies further parsing logic
    if (schemeLess.empty() || native_pathIsAbsolute(schemeLess))
        return schemeLess;

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
        if (schemeLess.find(native_pathSeparator) == 1)
            schemeLess = schemeLess.substr(2);
    }
    return from + native_pathSeparator + schemeLess;
}

#else // else for non-Windows systems

static const char native_pathSeparator = '/';

static bool native_pathIsAbsolute(const string& filePath)
{ return filePath.at(0) == '/'; }

static string native_resolvedFilePath(const string& filePath,
                                      const string& from)
{
    // strip file:// scheme if present
    auto schemeLess = removeFileUriScheme(filePath);

    // if empty or already absolute path, just return it
    // upfront empty check simplifies further parsing logic
    if (schemeLess.empty() || native_pathIsAbsolute(schemeLess))
        return schemeLess;

    // else make relative from the provided 'from' directory
    //
    // since we're prepending the 'from' directory, we can remove
    // any leading './' form our file path. this may just mean that
    // we pop it off to add it right back (when from == '.'), but this
    // keeps it consistent with other 'from' parent directories
    //
    const bool thisDirAtStart = (schemeLess.find(".") == 0);
    if (thisDirAtStart) {
        if (schemeLess.find(native_pathSeparator) == 1)
            schemeLess = schemeLess.substr(2);
    }
    return from + native_pathSeparator + schemeLess;
}

#endif // PBBAM_WIN_FILEPATHS

// see http://stackoverflow.com/questions/2869594/how-return-a-stdstring-from-cs-getcwd-function
string FileUtils::CurrentWorkingDirectory(void)
{
    const size_t chunkSize = 1024;
    const size_t maxNumChunks = 20;

    // stack-based buffer for 'normal' case
    char buffer[chunkSize];
    if (getcwd(buffer, sizeof(buffer)) != NULL)
        return string(buffer);

    // if error is not ERANGE, then it's not a problem of too-long name... something else happened
    if (errno != ERANGE)
        throw runtime_error("could not determine current working directory path");

    // long path - use heap, trying progressively longer buffers
    for (size_t chunks = 2; chunks < maxNumChunks; ++chunks) {
        unique_ptr<char> cwd(new char[chunkSize*chunks]);
        if (getcwd(cwd.get(), chunkSize*chunks) != NULL)
            return string(cwd.get());

        // if error is not ERANGE, then it's not a problem of too-long name... something else happened
        if (errno != ERANGE)
            throw runtime_error("could not determine current working directory path");
    }

    // crazy long path name
    throw runtime_error("could determine current working directory - extremely long path");
}

string FileUtils::DirectoryName(const string& file)
{
    const size_t found = file.rfind(Separator(), file.length());
    if (found != string::npos)
        return file.substr(0, found);
    return string(".");
}

bool FileUtils::Exists(const char* fn)
{
    struct stat buf;
    return (stat(fn, &buf) != -1);
}

chrono::system_clock::time_point FileUtils::LastModified(const char* fn)
{
    struct stat s;
    if (stat(fn, &s) != 0)
        throw runtime_error("could not get file timestamp");
    return chrono::system_clock::from_time_t(s.st_mtime);
}

string FileUtils::ResolvedFilePath(const string& filePath,
                                   const string& from)
{ return native_resolvedFilePath(filePath, from); }

constexpr char FileUtils::Separator(void)
{ return native_pathSeparator; }

off_t FileUtils::Size(const char* fn)
{
    struct stat s;
    if (stat(fn, &s) != 0)
        throw runtime_error("could not determine file size");
    return s.st_size;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio
