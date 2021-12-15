#ifndef PBBAM_TIMEUTILS_H
#define PBBAM_TIMEUTILS_H

#include <pbbam/Config.h>

#include <chrono>
#include <stdexcept>
#include <string>

#include <ctime>

namespace PacBio {
namespace BAM {

class TimeUtils
{
public:
    static std::string ToIso8601(const std::chrono::system_clock::time_point& tp)
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

    static std::string ToDataSetFormat(const std::chrono::system_clock::time_point& tp)
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
        if (ms.count() > 0) {
            result.append(std::to_string(ms.count()));
        }
        return result;
    }

    static std::chrono::system_clock::time_point CurrentTime()
    {
        return std::chrono::system_clock::now();
    }
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_TIMEUTILS_H
