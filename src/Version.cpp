// File Description
/// \file Version.cpp
/// \brief Implements the Version class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "Version.h"

#include <sstream>
#include <stdexcept>

#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {

const Version Version::Current = Version(3, 0, 7);
const Version Version::Minimum = Version(3, 0, 1);

// string must be "<major>.<minor>.<version>"
Version::Version(const std::string& v) : major_{0}, minor_{0}, revision_{0}
{
    // parse string
    try {
        const auto fields = Split(v, '.');
        const auto numFields = fields.size();
        if (numFields == 0) throw std::runtime_error{"Version: empty string"};
        major_ = std::stoi(fields.at(0));
        if (numFields > 1) {
            minor_ = std::stoi(fields.at(1));
            if (numFields > 2) revision_ = std::stoi(fields.at(2));
        }
    } catch (std::exception&) {
        throw std::runtime_error{"Version: could not parse: " + v};
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

}  // namespace BAM
}  // namespace PacBio
