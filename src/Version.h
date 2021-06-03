#ifndef PBBAM_VERSION_H
#define PBBAM_VERSION_H

#include <pbbam/Config.h>

#include <iosfwd>
#include <stdexcept>
#include <string>
#include <tuple>

namespace PacBio {
namespace BAM {

class Version
{
public:
    static const Version Current;
    static const Version Minimum;

    constexpr Version() = default;

    Version(int major, int minor, int revision);

    // string must be "<major>.<minor>.<version>"
    explicit Version(const std::string& v);

    Version(const Version&) = default;
    Version(Version&&) noexcept = default;
    Version& operator=(const Version&) = default;
    Version& operator=(Version&&) noexcept = default;
    ~Version() = default;

    bool operator==(const Version& other) const noexcept;
    bool operator!=(const Version& other) const noexcept;
    bool operator<(const Version& other) const noexcept;
    bool operator<=(const Version& other) const noexcept;
    bool operator>(const Version& other) const noexcept;
    bool operator>=(const Version& other) const noexcept;

    std::string ToString() const;
    explicit operator std::string() const;

    int Major() const;
    int Minor() const;
    int Revision() const;

    Version& Major(int major);
    Version& Minor(int minor);
    Version& Revision(int revision);

private:
    void Check() const;

    int major_ = 0;
    int minor_ = 0;
    int revision_ = 0;
};

std::ostream& operator<<(std::ostream& out, const Version& version);

inline Version::Version(int major, int minor, int revision)
    : major_{major}, minor_{minor}, revision_{revision}
{
    Check();
}

inline bool Version::operator==(const Version& other) const noexcept
{
    return std::tie(major_, minor_, revision_) ==
           std::tie(other.major_, other.minor_, other.revision_);
}

inline bool Version::operator!=(const Version& other) const noexcept { return !(*this == other); }

inline bool Version::operator<(const Version& other) const noexcept
{
    return std::tie(major_, minor_, revision_) <
           std::tie(other.major_, other.minor_, other.revision_);
}
inline bool Version::operator<=(const Version& other) const noexcept { return !(*this > other); }

inline bool Version::operator>(const Version& other) const noexcept { return other < *this; }

inline bool Version::operator>=(const Version& other) const noexcept { return !(*this < other); }

inline Version::operator std::string() const { return ToString(); }

inline void Version::Check() const
{
    if (major_ < 0 || minor_ < 0 || revision_ < 0) {
        throw std::runtime_error{"[pbbam] version string ERROR: cannot contain negative numbers"};
    }
}

inline int Version::Major() const { return major_; }

inline Version& Version::Major(int major)
{
    major_ = major;
    Check();
    return *this;
}

inline int Version::Minor() const { return minor_; }

inline Version& Version::Minor(int minor)
{
    minor_ = minor;
    Check();
    return *this;
}

inline int Version::Revision() const { return revision_; }

inline Version& Version::Revision(int revision)
{
    revision_ = revision;
    Check();
    return *this;
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_VERSION_H
