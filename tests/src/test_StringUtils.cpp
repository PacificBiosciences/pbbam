// Author: Derek Barnett

#include <gtest/gtest.h>

#include <pbbam/StringUtilities.h>

TEST(StringUtilsTest, BasicSplitWithDefaultDelim)
{
    const auto tokens = PacBio::BAM::Split("foo\tbar\tbaz");
    ASSERT_EQ(3, tokens.size());
    EXPECT_EQ("foo", tokens[0]);
    EXPECT_EQ("bar", tokens[1]);
    EXPECT_EQ("baz", tokens[2]);
}

TEST(StringUtilsTest, BasicSplitWithProvidedDelim)
{
    const auto tokens = PacBio::BAM::Split("foo:bar:baz", ':');
    ASSERT_EQ(3, tokens.size());
    EXPECT_EQ("foo", tokens[0]);
    EXPECT_EQ("bar", tokens[1]);
    EXPECT_EQ("baz", tokens[2]);
}

TEST(StringUtilsTest, SplitEmptyStringReturnsEmptyResult)
{
    const auto tokens = PacBio::BAM::Split("");
    EXPECT_TRUE(tokens.empty());
}

TEST(StringUtilsTest, SplitKeepsEmptyTokens)
{
    const auto tokens = PacBio::BAM::Split("foo\tbar\t\tbaz");
    ASSERT_EQ(4, tokens.size());
    EXPECT_EQ("foo", tokens[0]);
    EXPECT_EQ("bar", tokens[1]);
    EXPECT_EQ("", tokens[2]);
    EXPECT_EQ("baz", tokens[3]);
}

TEST(StringUtilsTest, RemoveWhitespaceNormal)
{
    const auto result =
        PacBio::BAM::RemoveAllWhitespace(" \f\r\v  Lorem ipsum     \tdolor sit\n\namet ");
    EXPECT_EQ("Loremipsumdolorsitamet", result);
}

TEST(StringUtilsTest, RemoveWhitespaceOnEmptyString)
{
    const auto result = PacBio::BAM::RemoveAllWhitespace("");
    EXPECT_TRUE(result.empty());
}
