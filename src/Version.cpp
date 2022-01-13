#include "PbbamInternalConfig.h"

#include "Version.h"

#include <pbbam/StringUtilities.h>

#include <ostream>
#include <sstream>
#include <stdexcept>

namespace PacBio {
namespace BAM {

const Version Version::Current = Version(5, 0, 0);
const Version Version::Minimum = Version(3, 0, 1);

// string must be "<major>.<minor>.<version>"
Version::Version(const std::string& v) : major_{0}, minor_{0}, revision_{0}
{
    // parse string
    try {
        const auto fields = Split(v, '.');
        const auto numFields = fields.size();
        if (numFields == 0) {
            throw std::runtime_error{"[pbbam] version string parsing ERROR: empty string"};
        }
        major_ = std::stoi(fields.at(0));
        if (numFields > 1) {
            minor_ = std::stoi(fields.at(1));
            if (numFields > 2) {
                revision_ = std::stoi(fields.at(2));
            }
        }
    } catch (std::exception& e) {
        std::ostringstream msg;
        msg << "[pbbam] version string parsing ERROR: failed to parse:\n"
            << "  version: " << v << '\n'
            << "  reason: " << e.what();
        throw std::runtime_error{msg.str()};
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

std::ostream& operator<<(std::ostream& out, const Version& version)
{
    out << version.ToString();
    return out;
}

}  // namespace BAM
}  // namespace PacBio
