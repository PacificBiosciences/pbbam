// Author: Derek Barnett

#ifndef TIMEUTILS_H
#define TIMEUTILS_H

#include <cassert>
#include <chrono>
#include <ctime>
#include <stdexcept>
#include <string>

namespace PacBio {
namespace BAM {
namespace internal {

inline std::string ToIso8601(const std::chrono::system_clock::time_point& tp)
{
    // get time info
    const time_t ttime_t = std::chrono::system_clock::to_time_t(tp);
    const std::chrono::system_clock::time_point tp_sec =
        std::chrono::system_clock::from_time_t(ttime_t);
    const std::chrono::milliseconds ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(tp - tp_sec);
    const std::tm* ttm =
        gmtime(&ttime_t);  // static obj, no free needed (may not be thread-safe though)

    // format output
    constexpr static const char date_time_format[] = "%FT%T";
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

inline std::string ToDataSetFormat(const std::chrono::system_clock::time_point& tp)
{
    // get time info
    const time_t ttime_t = std::chrono::system_clock::to_time_t(tp);
    const std::chrono::system_clock::time_point tp_sec =
        std::chrono::system_clock::from_time_t(ttime_t);
    const std::chrono::milliseconds ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(tp - tp_sec);
    const std::tm* ttm =
        gmtime(&ttime_t);  // static obj, no free needed (may not be thread-safe though)

    // format output
    constexpr static const char date_time_format[] = "%y%m%d_%H%M%S";
    char date_time_str[50];
    strftime(date_time_str, sizeof(date_time_str), date_time_format, ttm);
    std::string result(date_time_str);
    if (ms.count() > 0) result.append(std::to_string(ms.count()));
    return result;
}

inline std::chrono::system_clock::time_point CurrentTime()
{
    return std::chrono::system_clock::now();
}

}  // namespace PacBio
}  // namespace BAM
}  // namespace internal

#endif  // TIMEUTILS_H
