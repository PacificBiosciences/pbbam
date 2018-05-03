// File Description
/// \file Version.h
/// \brief Defines the Version class.
//
// Author: Derek Barnett

#ifndef PACBIOBAM_VERSION_H
#define PACBIOBAM_VERSION_H

#include <ostream>
#include <stdexcept>
#include <string>

namespace PacBio {
namespace BAM {
namespace internal {

class Version
{
public:
    static const Version Current;
    static const Version Minimum;

public:
    constexpr Version();

    Version(int major, int minor, int revision);

    // string must be "<major>.<minor>.<version>"
    Version(const std::string& v);

    Version(const Version&) = default;
    Version(Version&&) = default;
    Version& operator=(const Version&) = default;
    Version& operator=(Version&&) = default;
    ~Version() = default;

public:
    bool operator==(const Version& other) const;
    bool operator!=(const Version& other) const;
    bool operator<(const Version& other) const;
    bool operator<=(const Version& other) const;
    bool operator>(const Version& other) const;
    bool operator>=(const Version& other) const;

public:
    std::string ToString() const;
    operator std::string() const;

public:
    int Major() const;
    int Minor() const;
    int Revision() const;

public:
    Version& Major(int major);
    Version& Minor(int minor);
    Version& Revision(int revision);

private:
    int major_;
    int minor_;
    int revision_;

private:
    void Check() const;
};

inline std::ostream& operator<<(std::ostream& out, const Version& version)
{
    out << version.ToString();
    return out;
}

inline constexpr Version::Version() : major_{0}, minor_{0}, revision_{0} {}

inline Version::Version(int major, int minor, int revision)
    : major_{major}, minor_{minor}, revision_{revision}
{
    Check();
}

inline bool Version::operator==(const Version& other) const
{
    return major_ == other.major_ && minor_ == other.minor_ && revision_ == other.revision_;
}

inline bool Version::operator!=(const Version& other) const { return !(*this == other); }

inline bool Version::operator<(const Version& other) const
{
    // 2.* < 3.*
    if (major_ < other.major_) return true;

    // 3. ==  3.
    else if (major_ == other.major_) {

        // 3.1.* < 3.2.*
        if (minor_ < other.minor_) return true;

        // 3.2. == 3.2.
        else if (minor_ == other.minor_) {

            // 3.2.1 < 3.2.2
            if (revision_ < other.revision_) return true;
        }
    }

    // otherwise not less-than
    return false;
}
inline bool Version::operator<=(const Version& other) const { return !(*this > other); }

inline bool Version::operator>(const Version& other) const { return other < *this; }

inline bool Version::operator>=(const Version& other) const { return !(*this < other); }

inline Version::operator std::string() const { return ToString(); }

inline void Version::Check() const
{
    if (major_ < 0 || minor_ < 0 || revision_ < 0)
        throw std::runtime_error{"version cannot contain negative numbers"};
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

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // PACBIOBAM_VERSION_H
