// Author: Derek Barnett

#include <climits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/../../src/SequenceUtils.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;

TEST(SequenceUtilsTest, ComplementChar)
{
    // complement
    const char A = 'A';  // T
    const char B = 'B';  // V
    const char C = 'C';  // G
    const char D = 'D';  // H
    const char E = 'E';  // null
    const char F = 'F';  // null
    const char G = 'G';  // C
    const char H = 'H';  // D
    const char I = 'I';  // null
    const char J = 'J';  // null
    const char K = 'K';  // M
    const char L = 'L';  // null
    const char M = 'M';  // K
    const char N = 'N';  // N
    const char O = 'O';  // null
    const char P = 'P';  // null
    const char Q = 'Q';  // null
    const char R = 'R';  // Y
    const char S = 'S';  // S
    const char T = 'T';  // A
    const char U = 'U';  // A
    const char V = 'V';  // B
    const char W = 'W';  // W
    const char X = 'X';  // null
    const char Y = 'Y';  // R
    const char Z = 'Z';  // null

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
    std::string input1{"ATATATCCCGGCG"};
    const std::string rc1{"CGCCGGGATATAT"};

    ReverseComplement(input1);
    EXPECT_EQ(rc1, input1);
}
