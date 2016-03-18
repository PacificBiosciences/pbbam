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

// Author: Derek Barnett

#include "../src/Version.h"

#include <gtest/gtest.h>
#include <sstream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace tests {

static inline Version MakeVersion(int x, int y, int z)
{ return Version(x, y, z); }

} // namespace tests
} // namespace BAM
} // namespace PacBio

TEST(VersionTest, DefaultOk)
{
    Version v;
    EXPECT_EQ(0, v.Major());
    EXPECT_EQ(0, v.Minor());
    EXPECT_EQ(0, v.Revision());
}

TEST(VersionTest, CopyAndMoveOk)
{
    {   // copy ctor
        Version v1(3,1,1);
        EXPECT_EQ(3, v1.Major());
        EXPECT_EQ(1, v1.Minor());
        EXPECT_EQ(1, v1.Revision());

        Version v2(v1);
        EXPECT_EQ(3, v2.Major());
        EXPECT_EQ(1, v2.Minor());
        EXPECT_EQ(1, v2.Revision());
    }
    {   // copy assign
        Version v1(3,1,1);
        EXPECT_EQ(3, v1.Major());
        EXPECT_EQ(1, v1.Minor());
        EXPECT_EQ(1, v1.Revision());

        Version v2;
        v2 = v1;
        EXPECT_EQ(3, v2.Major());
        EXPECT_EQ(1, v2.Minor());
        EXPECT_EQ(1, v2.Revision());

    }
    {   // move ctor
        Version v(tests::MakeVersion(3,1,1));
        EXPECT_EQ(3, v.Major());
        EXPECT_EQ(1, v.Minor());
        EXPECT_EQ(1, v.Revision());

    }
    {   // move assign
        Version v1(3,1,1);
        EXPECT_EQ(3, v1.Major());
        EXPECT_EQ(1, v1.Minor());
        EXPECT_EQ(1, v1.Revision());

        Version v2;
        v2 = std::move(v1);
        EXPECT_EQ(3, v2.Major());
        EXPECT_EQ(1, v2.Minor());
        EXPECT_EQ(1, v2.Revision());
    }
}

TEST(VersionTest, FromIntsOk)
{
    {   // normal
        Version v(3,1,1);
        EXPECT_EQ(3, v.Major());
        EXPECT_EQ(1, v.Minor());
        EXPECT_EQ(1, v.Revision());
    }

    // negatives
    EXPECT_THROW(Version(-3, 1, 1), std::runtime_error);
}

TEST(VersionTest, FromStringOk)
{
    {   // normal
        Version v("3.1.1");
        EXPECT_EQ(3, v.Major());
        EXPECT_EQ(1, v.Minor());
        EXPECT_EQ(1, v.Revision());
    }

    // negatives
    EXPECT_THROW(Version("-3.1.1"), std::runtime_error);

    // non-numeric
    EXPECT_THROW(Version("foo.bar.baz"), std::runtime_error);

    // empty
    EXPECT_THROW(Version(""), std::runtime_error);
}

TEST(VersionTest, SettersOk)
{
    Version v(3,1,1);

    v.Major(4);

    EXPECT_EQ(4, v.Major());
    EXPECT_EQ(1, v.Minor());
    EXPECT_EQ(1, v.Revision());

    v.Minor(7);

    EXPECT_EQ(4, v.Major());
    EXPECT_EQ(7, v.Minor());
    EXPECT_EQ(1, v.Revision());

    v.Revision(23);

    EXPECT_EQ(4, v.Major());
    EXPECT_EQ(7, v.Minor());
    EXPECT_EQ(23, v.Revision());

    {   // invalid
        Version v1(3,1,1);
        Version v2(3,1,1);
        Version v3(3,1,1);
        EXPECT_THROW(v1.Major(-1),    std::runtime_error);
        EXPECT_THROW(v2.Minor(-1),    std::runtime_error);
        EXPECT_THROW(v3.Revision(-1), std::runtime_error);
    }
}

TEST(VersionTest, ComparisonsOk)
{
    const Version v0_0_0 = Version(0,0,0);
    const Version v0_0_4 = Version(0,0,4);
    const Version v0_1_0 = Version(0,1,0);
    const Version v0_1_4 = Version(0,1,4);
    const Version v3_0_0 = Version(3,0,0);
    const Version v3_0_4 = Version(3,0,4);
    const Version v3_1_0 = Version(3,1,0);
    const Version v3_1_4 = Version(3,1,4);
    const Version v3_1_5 = Version(3,1,5);

    // operator==
    EXPECT_TRUE(v0_0_0 == v0_0_0);
    EXPECT_TRUE(v3_0_0 == v3_0_0);
    EXPECT_TRUE(v0_1_0 == v0_1_0);
    EXPECT_TRUE(v0_0_4 == v0_0_4);
    EXPECT_TRUE(v3_1_0 == v3_1_0);
    EXPECT_TRUE(v3_1_4 == v3_1_4);

    EXPECT_FALSE(v3_1_4 == v0_0_0);
    EXPECT_FALSE(v3_1_4 == v3_0_0);
    EXPECT_FALSE(v3_1_4 == v0_1_0);
    EXPECT_FALSE(v3_1_4 == v0_0_4);
    EXPECT_FALSE(v3_1_4 == v3_1_0);
    EXPECT_FALSE(v3_1_4 == v3_1_5);

    // operator!=
    EXPECT_FALSE(v0_0_0 != v0_0_0);
    EXPECT_FALSE(v3_0_0 != v3_0_0);
    EXPECT_FALSE(v0_1_0 != v0_1_0);
    EXPECT_FALSE(v0_0_4 != v0_0_4);
    EXPECT_FALSE(v3_1_0 != v3_1_0);
    EXPECT_FALSE(v3_1_4 != v3_1_4);

    EXPECT_TRUE(v3_1_4 != v0_0_0);
    EXPECT_TRUE(v3_1_4 != v3_0_0);
    EXPECT_TRUE(v3_1_4 != v0_1_0);
    EXPECT_TRUE(v3_1_4 != v0_0_4);
    EXPECT_TRUE(v3_1_4 != v3_1_0);
    EXPECT_TRUE(v3_1_4 != v3_1_5);

    // operator<
    EXPECT_FALSE(v0_0_0 < v0_0_0);
    EXPECT_TRUE(v0_0_0 < v0_0_4);
    EXPECT_TRUE(v0_0_0 < v0_1_0);
    EXPECT_TRUE(v0_0_0 < v3_0_0);
    EXPECT_TRUE(v0_0_0 < v0_1_4);
    EXPECT_TRUE(v0_0_0 < v3_0_4);
    EXPECT_TRUE(v0_0_0 < v3_1_0);
    EXPECT_TRUE(v0_0_0 < v3_1_4);

    EXPECT_TRUE(v0_0_4 < v3_1_4);
    EXPECT_TRUE(v0_1_0 < v3_1_4);
    EXPECT_TRUE(v0_1_4 < v3_1_4);
    EXPECT_TRUE(v3_0_0 < v3_1_4);
    EXPECT_TRUE(v3_0_4 < v3_1_4);
    EXPECT_TRUE(v3_1_0 < v3_1_4);
    EXPECT_FALSE(v3_1_4 < v3_1_4);
    EXPECT_FALSE(v3_1_5 < v3_1_4);

    EXPECT_FALSE(v3_1_4 < v0_0_0);

    // operator<=
    EXPECT_TRUE(v0_0_0 <= v0_0_0);
    EXPECT_TRUE(v0_0_0 <= v0_0_4);
    EXPECT_TRUE(v0_0_0 <= v0_1_0);
    EXPECT_TRUE(v0_0_0 <= v3_0_0);
    EXPECT_TRUE(v0_0_0 <= v0_1_4);
    EXPECT_TRUE(v0_0_0 <= v3_0_4);
    EXPECT_TRUE(v0_0_0 <= v3_1_0);
    EXPECT_TRUE(v0_0_0 <= v3_1_4);

    EXPECT_TRUE(v0_0_4 <= v3_1_4);
    EXPECT_TRUE(v0_1_0 <= v3_1_4);
    EXPECT_TRUE(v0_1_4 <= v3_1_4);
    EXPECT_TRUE(v3_0_0 <= v3_1_4);
    EXPECT_TRUE(v3_0_4 <= v3_1_4);
    EXPECT_TRUE(v3_1_0 <= v3_1_4);
    EXPECT_TRUE(v3_1_4 <= v3_1_4);
    EXPECT_FALSE(v3_1_5 <= v3_1_4);

    EXPECT_FALSE(v3_1_4 <= v0_0_0);

    // operator>
    EXPECT_FALSE(v0_0_0 > v0_0_0);
    EXPECT_FALSE(v0_0_0 > v0_0_4);
    EXPECT_FALSE(v0_0_0 > v0_1_0);
    EXPECT_FALSE(v0_0_0 > v3_0_0);
    EXPECT_FALSE(v0_0_0 > v0_1_4);
    EXPECT_FALSE(v0_0_0 > v3_0_4);
    EXPECT_FALSE(v0_0_0 > v3_1_0);
    EXPECT_FALSE(v0_0_0 > v3_1_4);

    EXPECT_FALSE(v0_0_4 > v3_1_4);
    EXPECT_FALSE(v0_1_0 > v3_1_4);
    EXPECT_FALSE(v0_1_4 > v3_1_4);
    EXPECT_FALSE(v3_0_0 > v3_1_4);
    EXPECT_FALSE(v3_0_4 > v3_1_4);
    EXPECT_FALSE(v3_1_0 > v3_1_4);
    EXPECT_FALSE(v3_1_4 > v3_1_4);
    EXPECT_TRUE(v3_1_5 > v3_1_4);

    EXPECT_TRUE(v3_1_4 > v0_0_0);

    // operator>=
    EXPECT_TRUE(v0_0_0 >= v0_0_0);
    EXPECT_FALSE(v0_0_0 >= v0_0_4);
    EXPECT_FALSE(v0_0_0 >= v0_1_0);
    EXPECT_FALSE(v0_0_0 >= v3_0_0);
    EXPECT_FALSE(v0_0_0 >= v0_1_4);
    EXPECT_FALSE(v0_0_0 >= v3_0_4);
    EXPECT_FALSE(v0_0_0 >= v3_1_0);
    EXPECT_FALSE(v0_0_0 >= v3_1_4);

    EXPECT_FALSE(v0_0_4 >= v3_1_4);
    EXPECT_FALSE(v0_1_0 >= v3_1_4);
    EXPECT_FALSE(v0_1_4 >= v3_1_4);
    EXPECT_FALSE(v3_0_0 >= v3_1_4);
    EXPECT_FALSE(v3_0_4 >= v3_1_4);
    EXPECT_FALSE(v3_1_0 >= v3_1_4);
    EXPECT_TRUE(v3_1_4 >= v3_1_4);
    EXPECT_TRUE(v3_1_5 >= v3_1_4);

    EXPECT_TRUE(v3_1_4 >= v0_0_0);
}

TEST(VersionTest, ToStringOk)
{
    {
        Version v(0,0,0);
        EXPECT_EQ(string("0.0.0"), v.ToString());
    }
    {
        Version v(3,1,4);
        EXPECT_EQ(string("3.1.4"), v.ToString());
    }
    {
        Version v;
        v.Major(4);
        EXPECT_EQ(string("4.0.0"), v.ToString());
    }
    {
        const string s = "1.2.3";
        Version v(s);
        EXPECT_EQ(s, v.ToString());
    }
}

TEST(VersionTest, OutputStreamOk)
{
    Version v(3,1,4);
    Version v2(4,10,0);

    stringstream s;
    s << v << ", " << v2 << ", " << v << endl;

    EXPECT_EQ(string("3.1.4, 4.10.0, 3.1.4\n"), s.str());
}
