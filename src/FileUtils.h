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

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <sys/stat.h>

namespace PacBio {
namespace BAM {
namespace internal {

struct FileUtils
{
public:
    static bool Exists(const char* fn);
    static bool Exists(const std::string& fn);

    // throws if can't read
    static time_t LastModified(const char* fn);
    static time_t LastModified(const std::string& fn);

    // throws if can't read
    static off_t Size(const char* fn);
    static off_t Size(const std::string& fn);
};

inline bool FileUtils::Exists(const char* fn)
{ return Exists(std::string(fn)); }

inline bool FileUtils::Exists(const std::string& fn)
{
    std::ifstream stream(fn);
    return !stream.fail();
}

inline time_t FileUtils::LastModified(const char* fn)
{
    struct stat s;
    if (stat(fn, &s) != 0)
        throw std::exception();

#ifdef __DARWIN_64_BIT_INO_T
    return s.st_mtimespec.tv_sec; // 64-bit OSX has a modified stat struct
#else
    return s.st_mtime;            // all others?
#endif
}

inline time_t FileUtils::LastModified(const std::string& fn)
{ return LastModified(fn.c_str()); }

inline off_t FileUtils::Size(const char* fn)
{
    struct stat s;
    if (stat(fn, &s) != 0)
        throw std::exception();
    return s.st_size;
}

inline off_t FileUtils::Size(const std::string& fn)
{ return Size(fn.c_str()); }

//inline std::string FilenameExtension(const std::string& fn)
//{
//    const size_t lastDot = fn.find_last_of(".");
//    return (lastDot != std::string::npos ? fn.substr(lastDot+1) : std::string());
//}

////
//// -- examples --
////
//// input: /path/to/file.ext      result: file.ext
//// input: /path/to/file.ext.zip  result: file.ext.zip
//// input: file.ext               result: file.ext
////
//inline std::string FilenameFromPath(const std::string& fullPath)
//{
//    struct MatchesPathSeparator {
//        bool operator()(char c) const { return c == '/'; }
//    };

//    const auto lastSeparator = std::find_if(fullPath.rbegin(),
//                                            fullPath.rend(),
//                                            MatchesPathSeparator()).base();
//    return std::string(lastSeparator,fullPath.end());
//}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // FILEUTILS_H
