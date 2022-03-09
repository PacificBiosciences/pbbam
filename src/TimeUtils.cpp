#include <pbbam/Config.h>

#include <TimeUtils.h>

#include <chrono>
#include <ctime>
#include <mutex>
#include <string>

namespace PacBio {
namespace BAM {

std::string FormatTime(const time_t timeT, char const* dtFormat)
{
    char dateTimeStrBuf[50];  // NOLINT
    static std::mutex m{};
    {
        std::unique_lock<std::mutex> lk(m);
        const std::tm* ttm = gmtime(&timeT);  // NOLINT(concurrency-mt-unsafe)
        strftime(dateTimeStrBuf, sizeof(dateTimeStrBuf), dtFormat, ttm);
    }
    return {dateTimeStrBuf};
}

std::string TimeUtils::ToIso8601(const std::chrono::system_clock::time_point& tp)
{
    using namespace std::chrono;
    const time_t timeT = std::chrono::system_clock::to_time_t(tp);
    const system_clock::time_point tpSec = system_clock::from_time_t(timeT);
    const milliseconds ms = duration_cast<milliseconds>(tp - tpSec);

    auto result = FormatTime(timeT, "%FT%T");
    if (ms.count() > 0) {
        result.append(".");
        result.append(std::to_string(ms.count()));
    }
    result.append("Z");
    return result;
}

std::string TimeUtils::ToDataSetFormat(const std::chrono::system_clock::time_point& tp)
{
    using namespace std::chrono;
    const time_t timeT = std::chrono::system_clock::to_time_t(tp);
    const system_clock::time_point tpSec = system_clock::from_time_t(timeT);
    const milliseconds ms = duration_cast<milliseconds>(tp - tpSec);

    auto result = FormatTime(timeT, "%y%m%d_%H%M%S");
    if (ms.count() > 0) {
        result.append(std::to_string(ms.count()));
    }
    return result;
}

std::chrono::system_clock::time_point TimeUtils::CurrentTime()
{
    return std::chrono::system_clock::now();
}

}  // namespace BAM
}  // namespace PacBio
