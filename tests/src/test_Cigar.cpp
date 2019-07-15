// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include <pbbam/Cigar.h>

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;

TEST(CigarTest, TypeToCar)
{
    EXPECT_EQ('M', CigarOperation::TypeToChar(CigarOperationType::ALIGNMENT_MATCH) );
    EXPECT_EQ('I', CigarOperation::TypeToChar(CigarOperationType::INSERTION) );
    EXPECT_EQ('D', CigarOperation::TypeToChar(CigarOperationType::DELETION) );
    EXPECT_EQ('N', CigarOperation::TypeToChar(CigarOperationType::REFERENCE_SKIP) );
    EXPECT_EQ('S', CigarOperation::TypeToChar(CigarOperationType::SOFT_CLIP) );
    EXPECT_EQ('H', CigarOperation::TypeToChar(CigarOperationType::HARD_CLIP) );
    EXPECT_EQ('P', CigarOperation::TypeToChar(CigarOperationType::PADDING) );
    EXPECT_EQ('=', CigarOperation::TypeToChar(CigarOperationType::SEQUENCE_MATCH) );
    EXPECT_EQ('X', CigarOperation::TypeToChar(CigarOperationType::SEQUENCE_MISMATCH) );
}

TEST(CigarTest, CharToType)
{
    EXPECT_EQ(CigarOperationType::ALIGNMENT_MATCH,   CigarOperation::CharToType('M'));
    EXPECT_EQ(CigarOperationType::INSERTION,         CigarOperation::CharToType('I'));
    EXPECT_EQ(CigarOperationType::DELETION,          CigarOperation::CharToType('D'));
    EXPECT_EQ(CigarOperationType::REFERENCE_SKIP,    CigarOperation::CharToType('N'));
    EXPECT_EQ(CigarOperationType::SOFT_CLIP,         CigarOperation::CharToType('S'));
    EXPECT_EQ(CigarOperationType::HARD_CLIP,         CigarOperation::CharToType('H'));
    EXPECT_EQ(CigarOperationType::PADDING,           CigarOperation::CharToType('P'));
    EXPECT_EQ(CigarOperationType::SEQUENCE_MATCH,    CigarOperation::CharToType('='));
    EXPECT_EQ(CigarOperationType::SEQUENCE_MISMATCH, CigarOperation::CharToType('X'));
}

TEST(CigarTest, SetOperationYieldsCorrectType)
{
    CigarOperation c1; c1.Type(CigarOperationType::ALIGNMENT_MATCH);
    CigarOperation c2; c2.Type(CigarOperationType::INSERTION);
    CigarOperation c3; c3.Type(CigarOperationType::DELETION);
    CigarOperation c4; c4.Type(CigarOperationType::REFERENCE_SKIP);
    CigarOperation c5; c5.Type(CigarOperationType::SOFT_CLIP);
    CigarOperation c6; c6.Type(CigarOperationType::HARD_CLIP);
    CigarOperation c7; c7.Type(CigarOperationType::PADDING);
    CigarOperation c8; c8.Type(CigarOperationType::SEQUENCE_MATCH);
    CigarOperation c9; c9.Type(CigarOperationType::SEQUENCE_MISMATCH);

    EXPECT_EQ('M', c1.Char());
    EXPECT_EQ('I', c2.Char());
    EXPECT_EQ('D', c3.Char());
    EXPECT_EQ('N', c4.Char());
    EXPECT_EQ('S', c5.Char());
    EXPECT_EQ('H', c6.Char());
    EXPECT_EQ('P', c7.Char());
    EXPECT_EQ('=', c8.Char());
    EXPECT_EQ('X', c9.Char());
}

TEST(CigarTest, SetTypeYieldsCorrectOperation)
{
    CigarOperation c1; c1.Char('M');
    CigarOperation c2; c2.Char('I');
    CigarOperation c3; c3.Char('D');
    CigarOperation c4; c4.Char('N');
    CigarOperation c5; c5.Char('S');
    CigarOperation c6; c6.Char('H');
    CigarOperation c7; c7.Char('P');
    CigarOperation c8; c8.Char('=');
    CigarOperation c9; c9.Char('X');

    EXPECT_EQ(CigarOperationType::ALIGNMENT_MATCH,   c1.Type());
    EXPECT_EQ(CigarOperationType::INSERTION,         c2.Type());
    EXPECT_EQ(CigarOperationType::DELETION,          c3.Type());
    EXPECT_EQ(CigarOperationType::REFERENCE_SKIP,    c4.Type());
    EXPECT_EQ(CigarOperationType::SOFT_CLIP,         c5.Type());
    EXPECT_EQ(CigarOperationType::HARD_CLIP,         c6.Type());
    EXPECT_EQ(CigarOperationType::PADDING,           c7.Type());
    EXPECT_EQ(CigarOperationType::SEQUENCE_MATCH,    c8.Type());
    EXPECT_EQ(CigarOperationType::SEQUENCE_MISMATCH, c9.Type());
}

TEST(CigarStringTest, FromStdString_Empty)
{
    const std::string emptyCigar = "";
    Cigar cigar = Cigar::FromStdString(emptyCigar);
    EXPECT_TRUE(cigar.empty());
}

TEST(CigarStringTest, FromStdString_SingleOp)
{
    const std::string singleCigar = "100=";

    Cigar cigar = Cigar::FromStdString(singleCigar);
    ASSERT_TRUE(cigar.size() == 1);

    const CigarOperation& op = cigar.front();
    EXPECT_TRUE(op.Char()   == '=');
    EXPECT_TRUE(op.Length() == 100);
}

TEST(CigarStringTest, FromStdString_MultipleOps)
{
    const std::string multiCigar = "100=2D34I6=6X6=";

    Cigar cigar = Cigar::FromStdString(multiCigar);
    ASSERT_TRUE(cigar.size() == 6);

    CigarOperation op0 = cigar.at(0);
    CigarOperation op1 = cigar.at(1);
    CigarOperation op2 = cigar.at(2);
    CigarOperation op3 = cigar.at(3);
    CigarOperation op4 = cigar.at(4);
    CigarOperation op5 = cigar.at(5);

    EXPECT_TRUE(op0.Char()   == '=');
    EXPECT_TRUE(op0.Length() == 100);
    EXPECT_TRUE(op1.Char()   == 'D');
    EXPECT_TRUE(op1.Length() == 2);
    EXPECT_TRUE(op2.Char()   == 'I');
    EXPECT_TRUE(op2.Length() == 34);
    EXPECT_TRUE(op3.Char()   == '=');
    EXPECT_TRUE(op3.Length() == 6);
    EXPECT_TRUE(op4.Char()   == 'X');
    EXPECT_TRUE(op4.Length() == 6);
    EXPECT_TRUE(op5.Char()   == '=');
    EXPECT_TRUE(op5.Length() == 6);
}

TEST(CigarStringTest, ToStdString_Empty)
{
    const std::string empty;
    Cigar cigar;
    EXPECT_EQ(empty, cigar.ToStdString());
}

TEST(CigarStringTest, ToStdString_SingleOp)
{
    const std::string singleCigar = "100=";

    Cigar cigar;
    cigar.push_back( CigarOperation(CigarOperationType::SEQUENCE_MATCH, 100) );

    EXPECT_EQ(singleCigar, cigar.ToStdString());
}

TEST(CigarStringTest, ToStdString_MultipleOps)
{
    const std::string multiCigar = "100=2D34I6=6X6=";

    Cigar cigar;
    cigar.push_back(CigarOperation(CigarOperationType::SEQUENCE_MATCH,  100));
    cigar.push_back(CigarOperation(CigarOperationType::DELETION,          2));
    cigar.push_back(CigarOperation(CigarOperationType::INSERTION,        34));
    cigar.push_back(CigarOperation(CigarOperationType::SEQUENCE_MATCH,    6));
    cigar.push_back(CigarOperation(CigarOperationType::SEQUENCE_MISMATCH, 6));
    cigar.push_back(CigarOperation(CigarOperationType::SEQUENCE_MATCH,    6));

    EXPECT_EQ(multiCigar, cigar.ToStdString());
}

// clang-format on
