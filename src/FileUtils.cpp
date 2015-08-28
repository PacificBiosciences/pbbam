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
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

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
    char separator = '/';
#ifdef WIN32
    separator = '\\';
#endif
    const size_t found = file.rfind(separator, file.length());
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
{
    // strip scheme from beginning, if it exists
    string schemeLess = filePath;
    const size_t schemeFound = schemeLess.find("file://");
    if (schemeFound != string::npos) {
        if (schemeFound != 0)
            throw runtime_error("Malformed URI: scheme not at beginning");
        schemeLess = schemeLess.substr(7);
    }

    // if absolute path, just return it
    if (schemeLess.at(0) == '/')
        return schemeLess;
        
     // else make relative from the provided'from' directory (avoiding redundant "./"(s) at start)
    const bool thisDirAtStart = (schemeLess.find("./") == 0);
    if (thisDirAtStart)
        schemeLess = schemeLess.substr(2);
    return from + "/" + schemeLess;
}

off_t FileUtils::Size(const char* fn)
{
    struct stat s;
    if (stat(fn, &s) != 0)
        throw runtime_error("could not determine file size");
    return s.st_size;
}
