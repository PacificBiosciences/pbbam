#ifndef PBBAM_TAGCOLLECTION_H
#define PBBAM_TAGCOLLECTION_H

#include <pbbam/Config.h>

#include <pbbam/Tag.h>

#include <iosfwd>
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The TagCollection class represents a collection (or "dictionary") of
///        tags.
///
/// Tags are mapped to their tag name, a 2-character string.
///
class PBBAM_EXPORT TagCollection : public std::map<std::string, Tag>
{
public:
    /// \returns true if the collection contains a tag with \p name
    bool Contains(const std::string& name) const;

    /// \returns estimated number of bytes used by this tag collection
    ///
    /// \warning The actual usage is heavily implementation-dependent. A reasonable
    ///          estimate is provided here, but no guarantee can be made.
    ///
    int EstimatedBytesUsed() const noexcept;
};

///
/// Write tab-separated TagCollection to output stream. Tags are written in the form:
///
/// 'name=value'
///
std::ostream& operator<<(std::ostream& out, const TagCollection& tags);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_TAGCOLLECTION_H
