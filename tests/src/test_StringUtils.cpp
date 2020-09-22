// Author: Derek Barnett

#include <pbbam/StringUtilities.h>

#include <gtest/gtest.h>

TEST(BAM_StringUtils, can_split_using_default_delimiter)
{
    const auto tokens = PacBio::BAM::Split("foo\tbar\tbaz");
    ASSERT_EQ(3, tokens.size());
    EXPECT_EQ("foo", tokens[0]);
    EXPECT_EQ("bar", tokens[1]);
    EXPECT_EQ("baz", tokens[2]);
}

TEST(BAM_StringUtils, can_split_using_specified_delimiter)
{
    const auto tokens = PacBio::BAM::Split("foo:bar:baz", ':');
    ASSERT_EQ(3, tokens.size());
    EXPECT_EQ("foo", tokens[0]);
    EXPECT_EQ("bar", tokens[1]);
    EXPECT_EQ("baz", tokens[2]);
}

TEST(BAM_StringUtils, splitting_empty_string_returns_no_tokens)
{
    const auto tokens = PacBio::BAM::Split("");
    EXPECT_TRUE(tokens.empty());
}

TEST(BAM_StringUtils, splitting_with_consecutive_delimiters_keeps_empty_tokens)
{
    const auto tokens = PacBio::BAM::Split("foo\tbar\t\tbaz");
    ASSERT_EQ(4, tokens.size());
    EXPECT_EQ("foo", tokens[0]);
    EXPECT_EQ("bar", tokens[1]);
    EXPECT_EQ("", tokens[2]);
    EXPECT_EQ("baz", tokens[3]);
}

TEST(BAM_StringUtils, can_remove_whitespace_from_string)
{
    const auto result =
        PacBio::BAM::RemoveAllWhitespace(" \f\r\v  Lorem ipsum     \tdolor sit\n\namet ");
    EXPECT_EQ("Loremipsumdolorsitamet", result);
}

TEST(BAM_StringUtils, remove_whitespace_can_handle_empty_string)
{
    const auto result = PacBio::BAM::RemoveAllWhitespace("");
    EXPECT_TRUE(result.empty());
}
