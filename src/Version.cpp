// File Description
/// \file Version.cpp
/// \brief Implements the Version class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "Version.h"

#include <sstream>

#include "SequenceUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

const Version Version::Current = Version(3, 0, 5);
const Version Version::Minimum = Version(3, 0, 1);

// string must be "<major>.<minor>.<version>"
Version::Version(const std::string& v) : major_(0), minor_(0), revision_(0)
{
    // parse string
    try {
        const auto fields = internal::Split(v, '.');
        const auto numFields = fields.size();
        if (numFields == 0) throw std::runtime_error("invalid version number - empty string");
        major_ = std::stoi(fields.at(0));
        if (numFields > 1) {
            minor_ = std::stoi(fields.at(1));
            if (numFields > 2) revision_ = std::stoi(fields.at(2));
        }
    } catch (std::exception&) {
        auto msg = std::string{"invalid version number (" + v + "): failed to parse"};
        throw std::runtime_error(msg);
    }

    // ensure valid numbers
    Check();
}

std::string Version::ToString() const
{
    std::ostringstream s;
    s << major_ << '.' << minor_ << '.' << revision_;
    return s.str();
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
