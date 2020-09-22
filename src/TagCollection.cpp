#include "PbbamInternalConfig.h"

#include <pbbam/TagCollection.h>

#include <ostream>

namespace PacBio {
namespace BAM {

bool TagCollection::Contains(const std::string& name) const { return count(name) != 0; }

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
