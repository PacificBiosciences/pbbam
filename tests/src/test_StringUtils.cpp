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

#include <gtest/gtest.h>

#include <pbbam/StringUtilities.h>

TEST(StringUtilsTest, BasicSplitWithDefaultDelim)
{
    using PacBio::BAM::Split;

    const std::string test = "foo\tbar\tbaz";
    const auto tokens = Split(test);
    EXPECT_EQ(3, tokens.size());
    EXPECT_TRUE(tokens.at(0) == "foo");
    EXPECT_TRUE(tokens.at(1) == "bar");
    EXPECT_TRUE(tokens.at(2) == "baz");
}

TEST(StringUtilsTest, BasicSplitWithProvidedDelim)
{
    using PacBio::BAM::Split;

    const std::string test = "foo:bar:baz";
    const auto tokens = Split(test, ':');
    EXPECT_EQ(3, tokens.size());
    EXPECT_TRUE(tokens.at(0) == "foo");
    EXPECT_TRUE(tokens.at(1) == "bar");
    EXPECT_TRUE(tokens.at(2) == "baz");
}

TEST(StringUtilsTest, SplitEmptyStringReturnsEmptyResult)
{
    using PacBio::BAM::Split;

    const std::string test = "";
    const auto tokens = Split(test);
    EXPECT_TRUE(tokens.empty());
}

TEST(StringUtilsTest, SplitKeepsEmptyTokens)
{
    using PacBio::BAM::Split;

    const std::string test = "foo\tbar\t\tbaz";
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
        const std::string input = " \f\r\v  Lorem ipsum     \tdolor sit\n\namet ";
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
        const std::string input = "";
        const auto result = RemoveAllWhitespace(input);
        EXPECT_TRUE(result.empty());
    }
    {  // rvalue
        const auto result = RemoveAllWhitespace("");
        EXPECT_TRUE(result.empty());
    }
}
