// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
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
    constexpr Version(void);

    Version(int major, int minor, int revision);

    // string must be "<major>.<minor>.<version>"
    Version(const std::string& v);

    Version(const Version& other) = default;
    Version(Version&& other) = default;
    Version& operator=(const Version&) = default;
    Version& operator=(Version&&) = default;
    ~Version(void) = default;

public:
    bool operator==(const Version& other) const;
    bool operator!=(const Version& other) const;
    bool operator<(const Version& other) const;
    bool operator<=(const Version& other) const;
    bool operator>(const Version& other) const;
    bool operator>=(const Version& other) const;

public:
    std::string ToString(void) const;
    operator std::string(void) const;

public:
    int Major(void) const;
    int Minor(void) const;
    int Revision(void) const;

public:
    Version& Major(int major);
    Version& Minor(int minor);
    Version& Revision(int revision);

private:
    int major_;
    int minor_;
    int revision_;

private:
    void Check(void) const;
};

inline std::ostream& operator<<(std::ostream& out, const Version& version)
{
    out << version.ToString();
    return out;
}

inline constexpr Version::Version(void)
    : major_(0)
    , minor_(0)
    , revision_(0)
{ }

inline Version::Version(int major, int minor, int revision)
    : major_(major)
    , minor_(minor)
    , revision_(revision)
{ Check(); }

inline bool Version::operator==(const Version& other) const
{
    return major_ == other.major_ &&
           minor_ == other.minor_ &&
           revision_ == other.revision_;
}

inline bool Version::operator!=(const Version& other) const
{ return !(*this == other); }

inline bool Version::operator<(const Version& other) const
{
    // 2.* < 3.*
    if (major_ < other.major_)
        return true;

    // 3. ==  3.
    else if (major_ == other.major_) {

        // 3.1.* < 3.2.*
        if (minor_ < other.minor_)
            return true;

        // 3.2. == 3.2.
        else if (minor_ == other.minor_) {

            // 3.2.1 < 3.2.2
            if (revision_ < other.revision_)
                return true;
        }
    }

    // otherwise not less-than
    return false;
}
inline bool Version::operator<=(const Version& other) const
{ return !(*this > other); }

inline bool Version::operator>(const Version& other) const
{ return other < *this; }

inline bool Version::operator>=(const Version& other) const
{ return !(*this < other); }

inline Version::operator std::string(void) const
{ return ToString(); }

inline void Version::Check(void) const
{
    if (major_ < 0 || minor_ < 0 || revision_ < 0)
        throw std::runtime_error("version cannot contain negative numbers");
}

inline int Version::Major(void) const
{ return major_; }

inline Version& Version::Major(int major)
{
    major_ = major;
    Check();
    return *this;
}

inline int Version::Minor(void) const
{ return minor_; }

inline Version& Version::Minor(int minor)
{
    minor_ = minor;
    Check();
    return *this;
}

inline int Version::Revision(void) const
{ return revision_; }

inline Version& Version::Revision(int revision)
{
    revision_ = revision;
    Check();
    return *this;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // PACBIOBAM_VERSION_H
