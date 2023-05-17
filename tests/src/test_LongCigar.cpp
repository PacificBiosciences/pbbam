#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/StringUtilities.h>

#include "../../src/MemoryUtils.h"

#include "PbbamTestData.h"

using BamReader = PacBio::BAM::BamReader;
using BamRecord = PacBio::BAM::BamRecord;
using BamWriter = PacBio::BAM::BamWriter;
using Cigar = PacBio::Data::Cigar;
using CigarOp = PacBio::Data::CigarOperation;
using PacBio::Data::CigarOperationType;
using Tag = PacBio::BAM::Tag;

namespace LongCigarTests {

// BAM record in this file has its CIGAR data in the new "CG" tag
const std::string LongCigarBam = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/long-cigar-1.7.bam";

const std::string LongCigarOut =
    PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/long-cigar-generated.bam";

const std::size_t numOps = 72091;

BamRecord ReadLongCigarRecord(const std::string& fn)
{
    BamRecord b;
    BamReader reader{fn};
    const bool success = reader.GetNext(b);
    EXPECT_TRUE(success);
    return b;
}

}  // namespace LongCigarTests

TEST(BAM_LongCigar, can_read_long_cigar)
{
    const auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);

    EXPECT_EQ(LongCigarTests::numOps, b.CigarData().size());
    EXPECT_FALSE(b.Impl().HasTag("CG"));
}

TEST(BAM_LongCigar, can_edit_long_cigar)
{
    auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);
    b.Impl().CigarData(b.CigarData());

    EXPECT_EQ(LongCigarTests::numOps, b.CigarData().size());
    EXPECT_FALSE(b.Impl().HasTag("CG"));
}

TEST(BAM_LongCigar, can_write_long_cigar)
{
    {  // edit & write
        auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);
        b.Impl().CigarData(b.CigarData());

        EXPECT_EQ(LongCigarTests::numOps, b.CigarData().size());
        EXPECT_FALSE(b.Impl().HasTag("CG"));

        BamWriter writer{LongCigarTests::LongCigarOut, b.header_};
        writer.Write(b);
    }

    {  // read back in
        const auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarOut);

        EXPECT_EQ(LongCigarTests::numOps, b.CigarData().size());
        EXPECT_FALSE(b.Impl().HasTag("CG"));
    }
}
