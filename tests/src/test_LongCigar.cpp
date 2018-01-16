// Copyright (c) 2018, Pacific Biosciences of California, Inc.
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

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include "PbbamTestData.h"

#include <gtest/gtest.h>
#include <iostream>
#include <string>

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
