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

#ifdef PBBAM_TESTING
#define private public
#endif

#include <gtest/gtest.h>
#include <pbbam/../../src/SequenceUtils.h>
#include <pbbam/../../src/StringUtils.h>
#include <string>
#include <vector>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

TEST(StringUtilsTest, BasicSplit)
{
    const string test = "foo\tbar\tbaz";
    const vector<string> tokens = internal::Split(test, '\t');
    EXPECT_EQ(3, tokens.size());
    EXPECT_TRUE(tokens.at(0) == "foo");
    EXPECT_TRUE(tokens.at(1) == "bar");
    EXPECT_TRUE(tokens.at(2) == "baz");
}

TEST(StringUtilsTest, SplitKeepsEmptyTokens)
{
    const string test = "foo\tbar\t\tbaz";
    const vector<string> tokens = internal::Split(test, '\t');
    EXPECT_EQ(4, tokens.size());
    EXPECT_TRUE(tokens.at(0) == "foo");
    EXPECT_TRUE(tokens.at(1) == "bar");
    EXPECT_TRUE(tokens.at(2) == "");
    EXPECT_TRUE(tokens.at(3) == "baz");
}

/*
TEST(StringUtilsTest, Int2String)
{
    const int zero     = 0;
    const int one      = 1;
    const int negative = -42;
    const int large    = INT32_MAX;

    const string zeroString     = "0";
    const string oneString      = "1";
    const string negativeString = "-42";
    const string largeString    = "2147483647";

    EXPECT_NO_THROW(
    {
        EXPECT_EQ(zeroString,     internal::Int2String(zero));
        EXPECT_EQ(oneString,      internal::Int2String(one));
        EXPECT_EQ(negativeString, internal::Int2String(negative));
        EXPECT_EQ(largeString,    internal::Int2String(large));
    });
}


TEST(StringUtilsTest, String2Int)
{
    const string zeroString     = "0";
    const string oneString      = "1";
    const string negativeString = "-42";
    const string maxString      = "2147483647";
    const string minString      = "-2147483648";

    const string overflowString   = "2147483648";
    const string underflowString  = "-2147483649";
    const string floatNumString   = "3.14";
    const string notANumberString = "foo";

    const int zero     = 0;
    const int one      = 1;
    const int negative = -42;
    const int min      = INT32_MIN;
    const int max      = INT32_MAX;

    EXPECT_NO_THROW(
    {
        EXPECT_EQ(zero,     internal::String2Int(zeroString));
        EXPECT_EQ(one,      internal::String2Int(oneString));
        EXPECT_EQ(negative, internal::String2Int(negativeString));
        EXPECT_EQ(min,      internal::String2Int(minString));
        EXPECT_EQ(max,      internal::String2Int(maxString));
    });

    EXPECT_THROW(internal::String2Int(overflowString),   std::exception);
    EXPECT_THROW(internal::String2Int(underflowString),  std::exception);
    EXPECT_THROW(internal::String2Int(floatNumString),   std::exception);
    EXPECT_THROW(internal::String2Int(notANumberString), std::exception);
}
*/

TEST(SequenceUtilsTest, ComplementChar)
{
                        // complement
    const char A = 'A'; // T
    const char B = 'B'; // V
    const char C = 'C'; // G
    const char D = 'D'; // H
    const char E = 'E'; // null
    const char F = 'F'; // null
    const char G = 'G'; // C
    const char H = 'H'; // D
    const char I = 'I'; // null
    const char J = 'J'; // null
    const char K = 'K'; // M
    const char L = 'L'; // null
    const char M = 'M'; // K
    const char N = 'N'; // N
    const char O = 'O'; // null
    const char P = 'P'; // null
    const char Q = 'Q'; // null
    const char R = 'R'; // Y
    const char S = 'S'; // S
    const char T = 'T'; // A
    const char U = 'U'; // A
    const char V = 'V'; // B
    const char W = 'W'; // W
    const char X = 'X'; // null
    const char Y = 'Y'; // R
    const char Z = 'Z'; // null

    EXPECT_EQ(T, Complement(A));
    EXPECT_EQ(V, Complement(B));
    EXPECT_EQ(G, Complement(C));
    EXPECT_EQ(H, Complement(D));
    EXPECT_EQ(0, Complement(E));
    EXPECT_EQ(0, Complement(F));
    EXPECT_EQ(C, Complement(G));
    EXPECT_EQ(D, Complement(H));
    EXPECT_EQ(0, Complement(I));
    EXPECT_EQ(0, Complement(J));
    EXPECT_EQ(M, Complement(K));
    EXPECT_EQ(0, Complement(L));
    EXPECT_EQ(K, Complement(M));
    EXPECT_EQ(N, Complement(N));
    EXPECT_EQ(0, Complement(O));
    EXPECT_EQ(0, Complement(P));
    EXPECT_EQ(0, Complement(Q));
    EXPECT_EQ(Y, Complement(R));
    EXPECT_EQ(S, Complement(S));
    EXPECT_EQ(A, Complement(T));
    EXPECT_EQ(A, Complement(U));
    EXPECT_EQ(B, Complement(V));
    EXPECT_EQ(W, Complement(W));
    EXPECT_EQ(0, Complement(X));
    EXPECT_EQ(R, Complement(Y));
    EXPECT_EQ(0, Complement(Z));
}

TEST(SequenceUtilsTest, ReverseComplement)
{
    string input1 = "ATATATCCCGGCG";
    const string rc1 = "CGCCGGGATATAT";

    ReverseComplement(input1);
    EXPECT_EQ(rc1, input1);
}
