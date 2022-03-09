#ifndef PBBAM_TIMEUTILS_H
#define PBBAM_TIMEUTILS_H

#include <pbbam/Config.h>

#include <chrono>
#include <string>

namespace PacBio {
namespace BAM {

class TimeUtils
{
public:
    static std::string ToIso8601(const std::chrono::system_clock::time_point& tp);

    static std::string ToDataSetFormat(const std::chrono::system_clock::time_point& tp);

    static std::chrono::system_clock::time_point CurrentTime();
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_TIMEUTILS_H
