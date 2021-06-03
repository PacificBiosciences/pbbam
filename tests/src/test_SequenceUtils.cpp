#include <pbbam/../../src/SequenceUtils.h>

#include <climits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

TEST(BAM_SequenceUtils, can_complement_single_char)
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

    EXPECT_EQ(T, PacBio::BAM::Complement(A));
    EXPECT_EQ(V, PacBio::BAM::Complement(B));
    EXPECT_EQ(G, PacBio::BAM::Complement(C));
    EXPECT_EQ(H, PacBio::BAM::Complement(D));
    EXPECT_EQ(0, PacBio::BAM::Complement(E));
    EXPECT_EQ(0, PacBio::BAM::Complement(F));
    EXPECT_EQ(C, PacBio::BAM::Complement(G));
    EXPECT_EQ(D, PacBio::BAM::Complement(H));
    EXPECT_EQ(0, PacBio::BAM::Complement(I));
    EXPECT_EQ(0, PacBio::BAM::Complement(J));
    EXPECT_EQ(M, PacBio::BAM::Complement(K));
    EXPECT_EQ(0, PacBio::BAM::Complement(L));
    EXPECT_EQ(K, PacBio::BAM::Complement(M));
    EXPECT_EQ(N, PacBio::BAM::Complement(N));
    EXPECT_EQ(0, PacBio::BAM::Complement(O));
    EXPECT_EQ(0, PacBio::BAM::Complement(P));
    EXPECT_EQ(0, PacBio::BAM::Complement(Q));
    EXPECT_EQ(Y, PacBio::BAM::Complement(R));
    EXPECT_EQ(S, PacBio::BAM::Complement(S));
    EXPECT_EQ(A, PacBio::BAM::Complement(T));
    EXPECT_EQ(A, PacBio::BAM::Complement(U));
    EXPECT_EQ(B, PacBio::BAM::Complement(V));
    EXPECT_EQ(W, PacBio::BAM::Complement(W));
    EXPECT_EQ(0, PacBio::BAM::Complement(X));
    EXPECT_EQ(R, PacBio::BAM::Complement(Y));
    EXPECT_EQ(0, PacBio::BAM::Complement(Z));
}

TEST(BAM_SequenceUtils, can_reverse_complement_string)
{
    std::string seq{"ATATATCCCGGCG"};
    const std::string revComp{"CGCCGGGATATAT"};

    PacBio::BAM::ReverseComplement(seq);
    EXPECT_EQ(revComp, seq);
}
