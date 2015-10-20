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

#ifndef TIMEUTILS_H
#define TIMEUTILS_H

#include <chrono>
#include <stdexcept>
#include <string>
#include <cassert>
#include <ctime>

namespace PacBio {
namespace BAM {
namespace internal {

inline
std::string ToIso8601(const std::chrono::system_clock::time_point& tp)
{
    // get time info
    const time_t ttime_t = std::chrono::system_clock::to_time_t(tp);
    const std::chrono::system_clock::time_point tp_sec = std::chrono::system_clock::from_time_t(ttime_t);
    const std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp - tp_sec);
    const std::tm* ttm = gmtime(&ttime_t);  // static obj, no free needed (may not be thread-safe though)

    // format output
    char date_time_format[] = "%FT%T";
    char date_time_str[50];
    strftime(date_time_str, sizeof(date_time_str), date_time_format, ttm);
    std::string result(date_time_str);
    if (ms.count() > 0) {
        result.append(".");
        result.append(std::to_string(ms.count()));
    }
    result.append("Z");
    return result;
}

inline
std::string ToDataSetFormat(const std::chrono::system_clock::time_point& tp)
{
    // get time info
    const time_t ttime_t = std::chrono::system_clock::to_time_t(tp);
    const std::chrono::system_clock::time_point tp_sec = std::chrono::system_clock::from_time_t(ttime_t);
    const std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp - tp_sec);
    const std::tm* ttm = gmtime(&ttime_t);  // static obj, no free needed (may not be thread-safe though)

    // format output
    char date_time_format[] = "%y%m%d_%H%M%S";
    char date_time_str[50];
    strftime(date_time_str, sizeof(date_time_str), date_time_format, ttm);
    std::string result(date_time_str);
    if (ms.count() > 0)
        result.append(std::to_string(ms.count()));
    return result;
}

inline
std::chrono::system_clock::time_point CurrentTime(void)
{ return std::chrono::system_clock::now(); }

} // namespace PacBio
} // namespace BAM
} // namespace internal

#endif // TIMEUTILS_H
