// Author: Derek Barnett

#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "../src/Version.h"

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;

namespace VersionTests {

static inline Version MakeVersion(int x, int y, int z) { return Version(x, y, z); }

}  // namespace VersionTests

TEST(VersionTest, DefaultOk)
{
    Version v;
    EXPECT_EQ(0, v.Major());
    EXPECT_EQ(0, v.Minor());
    EXPECT_EQ(0, v.Revision());
}

TEST(VersionTest, CopyAndMoveOk)
{
    {  // copy ctor
        Version v1(3, 1, 1);
        EXPECT_EQ(3, v1.Major());
        EXPECT_EQ(1, v1.Minor());
        EXPECT_EQ(1, v1.Revision());

        Version v2(v1);
        EXPECT_EQ(3, v2.Major());
        EXPECT_EQ(1, v2.Minor());
        EXPECT_EQ(1, v2.Revision());
    }
    {  // copy assign
        Version v1(3, 1, 1);
        EXPECT_EQ(3, v1.Major());
        EXPECT_EQ(1, v1.Minor());
        EXPECT_EQ(1, v1.Revision());

        Version v2;
        v2 = v1;
        EXPECT_EQ(3, v2.Major());
        EXPECT_EQ(1, v2.Minor());
        EXPECT_EQ(1, v2.Revision());
    }
    {  // move ctor
        Version v(VersionTests::MakeVersion(3, 1, 1));
        EXPECT_EQ(3, v.Major());
        EXPECT_EQ(1, v.Minor());
        EXPECT_EQ(1, v.Revision());
    }
    {  // move assign
        Version v1(3, 1, 1);
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
    {  // normal
        Version v(3, 1, 1);
        EXPECT_EQ(3, v.Major());
        EXPECT_EQ(1, v.Minor());
        EXPECT_EQ(1, v.Revision());
    }

    // negatives
    EXPECT_THROW(Version(-3, 1, 1), std::runtime_error);
}

TEST(VersionTest, FromStringOk)
{
    {  // normal
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
    Version v(3, 1, 1);

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

    {  // invalid
        Version v1(3, 1, 1);
        Version v2(3, 1, 1);
        Version v3(3, 1, 1);
        EXPECT_THROW(v1.Major(-1), std::runtime_error);
        EXPECT_THROW(v2.Minor(-1), std::runtime_error);
        EXPECT_THROW(v3.Revision(-1), std::runtime_error);
    }
}

TEST(VersionTest, ComparisonsOk)
{
    const Version v0_0_0 = Version(0, 0, 0);
    const Version v0_0_4 = Version(0, 0, 4);
    const Version v0_1_0 = Version(0, 1, 0);
    const Version v0_1_4 = Version(0, 1, 4);
    const Version v3_0_0 = Version(3, 0, 0);
    const Version v3_0_4 = Version(3, 0, 4);
    const Version v3_1_0 = Version(3, 1, 0);
    const Version v3_1_4 = Version(3, 1, 4);
    const Version v3_1_5 = Version(3, 1, 5);

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
        Version v(0, 0, 0);
        EXPECT_EQ(std::string("0.0.0"), v.ToString());
    }
    {
        Version v(3, 1, 4);
        EXPECT_EQ(std::string("3.1.4"), v.ToString());
    }
    {
        Version v;
        v.Major(4);
        EXPECT_EQ(std::string("4.0.0"), v.ToString());
    }
    {
        const std::string s = "1.2.3";
        Version v(s);
        EXPECT_EQ(s, v.ToString());
    }
}

TEST(VersionTest, OutputStreamOk)
{
    Version v(3, 1, 4);
    Version v2(4, 10, 0);

    std::ostringstream s;
    s << v << ", " << v2 << ", " << v << std::endl;

    EXPECT_EQ(std::string("3.1.4, 4.10.0, 3.1.4\n"), s.str());
}
