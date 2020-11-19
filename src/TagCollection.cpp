#include "PbbamInternalConfig.h"

#include <pbbam/TagCollection.h>

#include <ostream>

namespace PacBio {
namespace BAM {

bool TagCollection::Contains(const std::string& name) const { return count(name) != 0; }

int TagCollection::EstimatedBytesUsed() const noexcept
{
    int result = sizeof(std::map<std::string, Tag>);
    for (const auto& tag : *this) {
        // estimated red-black tree node overhead: 3 ptrs + color enum
        result += (3 * sizeof(std::ptrdiff_t)) + sizeof(int);
        // node payload
        result += sizeof(std::string);  // 2-char tag names fit into SSO
        result += tag.second.EstimatedBytesUsed();
    }
    return result;
}

std::ostream& operator<<(std::ostream& out, const TagCollection& tags)
{
    bool first = true;
    for (const auto& tag : tags) {
        if (!first)
            out << '\t';
        else
            first = false;
        out << tag.first << '=' << tag.second;
    }

    return out;
}

}  // namespace BAM
}  // namespace PacBio
