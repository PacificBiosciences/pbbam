#include "../src/Version.h"

#include <sstream>
#include <string>

#include <gtest/gtest.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace VersionTests {

Version MakeVersion(int x, int y, int z) { return Version{x, y, z}; }

}  // namespace VersionTests

TEST(BAM_Version, default_all_zeros)
{
    Version v;
    EXPECT_EQ(0, v.Major());
    EXPECT_EQ(0, v.Minor());
    EXPECT_EQ(0, v.Revision());
}

TEST(BAM_Version, can_create_from_integers)
{
    Version v{3, 1, 1};
    EXPECT_EQ(3, v.Major());
    EXPECT_EQ(1, v.Minor());
    EXPECT_EQ(1, v.Revision());
}

TEST(BAM_Version, throws_on_negative_integers)
{
    EXPECT_THROW(Version(-3, 1, 1), std::runtime_error);
}

TEST(BAM_Version, can_create_from_string)
{
    Version v{"3.1.1"};
    EXPECT_EQ(3, v.Major());
    EXPECT_EQ(1, v.Minor());
    EXPECT_EQ(1, v.Revision());
}

TEST(BAM_Version, throws_on_negative_in_string)
{
    EXPECT_THROW(Version("-3.1.1"), std::runtime_error);
}

TEST(BAM_Version, throws_on_non_numeric_string)
{
    EXPECT_THROW(Version("foo.bar.baz"), std::runtime_error);
}

TEST(BAM_Version, throws_on_empty_string) { EXPECT_THROW(Version(""), std::runtime_error); }

TEST(BAM_Version, can_be_compared)
{
    const Version v0_0_0{0, 0, 0};
    const Version v0_0_4{0, 0, 4};
    const Version v0_1_0{0, 1, 0};
    const Version v0_1_4{0, 1, 4};
    const Version v3_0_0{3, 0, 0};
    const Version v3_0_4{3, 0, 4};
    const Version v3_1_0{3, 1, 0};
    const Version v3_1_4{3, 1, 4};
    const Version v3_1_5{3, 1, 5};

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

TEST(BAM_Version, can_be_converted_to_string)
{
    {
        const Version v{0, 0, 0};
        EXPECT_EQ("0.0.0", v.ToString());
    }
    {
        const Version v{3, 1, 4};
        EXPECT_EQ("3.1.4", v.ToString());
    }
    {
        Version v;
        v.Major(4);
        EXPECT_EQ("4.0.0", v.ToString());
    }
    {
        const std::string s{"1.2.3"};
        const Version v{s};
        EXPECT_EQ(s, v.ToString());
    }
}

TEST(BAM_Version, can_write_to_ostream)
{
    const Version v{3, 1, 4};
    const Version v2{4, 10, 0};

    std::ostringstream s;
    s << v << ", " << v2 << ", " << v << '\n';

    EXPECT_EQ("3.1.4, 4.10.0, 3.1.4\n", s.str());
}
