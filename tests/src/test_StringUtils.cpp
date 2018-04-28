// Author: Derek Barnett

#include <gtest/gtest.h>

#include <pbbam/StringUtilities.h>

TEST(StringUtilsTest, BasicSplitWithDefaultDelim)
{
    using PacBio::BAM::Split;

    const std::string test{"foo\tbar\tbaz"};
    const auto tokens = Split(test);
    EXPECT_EQ(3, tokens.size());
    EXPECT_TRUE(tokens.at(0) == "foo");
    EXPECT_TRUE(tokens.at(1) == "bar");
    EXPECT_TRUE(tokens.at(2) == "baz");
}

TEST(StringUtilsTest, BasicSplitWithProvidedDelim)
{
    using PacBio::BAM::Split;

    const std::string test{"foo:bar:baz"};
    const auto tokens = Split(test, ':');
    EXPECT_EQ(3, tokens.size());
    EXPECT_TRUE(tokens.at(0) == "foo");
    EXPECT_TRUE(tokens.at(1) == "bar");
    EXPECT_TRUE(tokens.at(2) == "baz");
}

TEST(StringUtilsTest, SplitEmptyStringReturnsEmptyResult)
{
    using PacBio::BAM::Split;

    const std::string test;
    const auto tokens = Split(test);
    EXPECT_TRUE(tokens.empty());
}

TEST(StringUtilsTest, SplitKeepsEmptyTokens)
{
    using PacBio::BAM::Split;

    const std::string test{"foo\tbar\t\tbaz"};
    const auto tokens = Split(test);
    EXPECT_EQ(4, tokens.size());
    EXPECT_TRUE(tokens.at(0) == "foo");
    EXPECT_TRUE(tokens.at(1) == "bar");
    EXPECT_TRUE(tokens.at(2) == "");
    EXPECT_TRUE(tokens.at(3) == "baz");
}

TEST(StringUtilsTest, RemoveWhitespaceNormal)
{
    using PacBio::BAM::RemoveAllWhitespace;

    {  // lvalue
        const std::string input{" \f\r\v  Lorem ipsum     \tdolor sit\n\namet "};
        const auto result = RemoveAllWhitespace(input);
        EXPECT_EQ("Loremipsumdolorsitamet", result);
    }
    {  // rvalue
        const auto result = RemoveAllWhitespace(" \f\r\v  Lorem ipsum     \tdolor sit\n\namet ");
        EXPECT_EQ("Loremipsumdolorsitamet", result);
    }
}

TEST(StringUtilsTest, RemoveWhitespaceOnEmptyString)
{
    using PacBio::BAM::RemoveAllWhitespace;

    {  // lvalue
        const std::string input;
        const auto result = RemoveAllWhitespace(input);
        EXPECT_TRUE(result.empty());
    }
    {  // rvalue
        const auto result = RemoveAllWhitespace("");
        EXPECT_TRUE(result.empty());
    }
}
