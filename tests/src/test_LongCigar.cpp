// Author: Derek Barnett

#include <iostream>
#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>

using BamReader = PacBio::BAM::BamReader;
using BamRecord = PacBio::BAM::BamRecord;
using BamWriter = PacBio::BAM::BamWriter;
using Cigar = PacBio::BAM::Cigar;
using CigarOp = PacBio::BAM::CigarOperation;
using PacBio::BAM::CigarOperationType;
using Tag = PacBio::BAM::Tag;

namespace LongCigarTests {

// BAM record in this file has its CIGAR data in the new "CG" tag
static const std::string LongCigarBam = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/long-cigar.bam";

static const std::string LongCigarOut =
    PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/long-cigar-generated.bam";

static const size_t numOps = 66000;

static BamRecord ReadLongCigarRecord(const std::string& fn)
{
    BamRecord b;
    BamReader reader{fn};
    const bool success = reader.GetNext(b);
    EXPECT_TRUE(success);
    return b;
}

static void SetLongCigar(BamRecord* b)
{
    Cigar cigar;
    cigar.resize(numOps);
    for (size_t i = 0; i < LongCigarTests::numOps; ++i) {
        const CigarOperationType type =
            (i % 2 == 0 ? CigarOperationType::SEQUENCE_MATCH : CigarOperationType::INSERTION);
        cigar.at(i) = CigarOp(type, 1);
    }
    b->Impl().CigarData(cigar);
}

static void CheckLongCigar(const Cigar& cigar)
{
    ASSERT_EQ(numOps, cigar.size());

    for (size_t i = 0; i < numOps; ++i) {
        const CigarOp& op = cigar.at(i);
        EXPECT_EQ(1, op.Length());

        const CigarOperationType expectedType =
            (i % 2 == 0 ? CigarOperationType::SEQUENCE_MATCH : CigarOperationType::INSERTION);
        EXPECT_EQ(expectedType, op.Type());
    }
}

static void CheckLongCigarTag(const Tag& cigarTag)
{
    ASSERT_TRUE(cigarTag.IsUInt32Array());
    const auto tagArray = cigarTag.ToUInt32Array();
    ASSERT_EQ(numOps, tagArray.size());

    for (size_t i = 0; i < numOps; ++i) {
        const auto op = tagArray.at(i);
        const auto expectedLength = 1;
        const auto expectedType = (i % 2 == 0 ? BAM_CEQUAL : BAM_CINS);

        EXPECT_EQ(expectedType, bam_cigar_op(op));
        EXPECT_EQ(expectedLength, bam_cigar_oplen(op));
    }
}

}  // namespace LongCigarTests

TEST(LongCigarTest, ReadAndFetchLongCigar)
{
    const auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);

    // public API
    const auto cigar = b.CigarData();
    EXPECT_EQ(66000, cigar.size());

    // TODO: come back & check raw data once we have 'private access wrapper'
    //       but we're looking good
}

TEST(LongCigarTest, EditLongCigar)
{
    SCOPED_TRACE("EditLongCigar");

    auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);
    LongCigarTests::SetLongCigar(&b);

    const auto recordCigar = b.CigarData();
    const auto cigarTag = b.Impl().TagValue("CG");
    LongCigarTests::CheckLongCigar(recordCigar);
    LongCigarTests::CheckLongCigarTag(cigarTag);
}

TEST(LongCigarTest, WriteLongCigar)
{
    SCOPED_TRACE("WriteLongCigar");

    {  // write record with our custom long CIGAR
        auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);
        LongCigarTests::SetLongCigar(&b);
        BamWriter writer{LongCigarTests::LongCigarOut, b.header_};
        writer.Write(b);
    }
    {  // read back in to check
        auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarOut);
        const auto recordCigar = b.CigarData();
        const auto cigarTag = b.Impl().TagValue("CG");
        LongCigarTests::CheckLongCigar(recordCigar);
        LongCigarTests::CheckLongCigarTag(cigarTag);
    }
}
