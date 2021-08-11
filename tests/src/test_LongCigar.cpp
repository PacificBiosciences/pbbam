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
using Cigar = PacBio::BAM::Cigar;
using CigarOp = PacBio::BAM::CigarOperation;
using PacBio::BAM::CigarOperationType;
using Tag = PacBio::BAM::Tag;

namespace LongCigarTests {

bool DoesHtslibSupportLongCigar()
{
    const std::string htsVersion = hts_version();

    // remove any "-<blah>" for non-release versions
    const auto versionBase = PacBio::BAM::Split(htsVersion, '-');
    if (versionBase.empty()) {
        throw std::runtime_error{"invalid htslib version format: " + htsVersion};
    }

    // grab major/minor version numbers
    const auto versionParts = PacBio::BAM::Split(versionBase[0], '.');
    if (versionParts.size() < 2) {
        throw std::runtime_error{"invalid htslib version format: " + htsVersion};
    }

    // check against v1.7
    const int versionMajor = std::stoi(versionParts[0]);
    const int versionMinor = std::stoi(versionParts[1]);
    static constexpr const int v17_major = 1;
    static constexpr const int v17_minor = 7;
    return std::tie(versionMajor, versionMinor) >= std::tie(v17_major, v17_minor);
}

const bool has_native_long_cigar_support = DoesHtslibSupportLongCigar();

// BAM record in this file has its CIGAR data in the new "CG" tag
const std::string LongCigarBam = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/long-cigar-1.7.bam";

const std::string LongCigarOut =
    PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/long-cigar-generated.bam";

const size_t numOps = 72091;

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
    if (LongCigarTests::has_native_long_cigar_support) {
        EXPECT_FALSE(b.Impl().HasTag("CG"));
    } else {
        EXPECT_TRUE(b.Impl().HasTag("CG"));
    }
}

TEST(BAM_LongCigar, can_edit_long_cigar)
{
    auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);
    b.Impl().CigarData(b.CigarData());

    EXPECT_EQ(LongCigarTests::numOps, b.CigarData().size());
    if (LongCigarTests::has_native_long_cigar_support) {
        EXPECT_FALSE(b.Impl().HasTag("CG"));
    } else {
        EXPECT_TRUE(b.Impl().HasTag("CG"));
    }
}

TEST(BAM_LongCigar, can_write_long_cigar)
{
    {  // edit & write
        auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarBam);
        b.Impl().CigarData(b.CigarData());

        EXPECT_EQ(LongCigarTests::numOps, b.CigarData().size());
        if (LongCigarTests::has_native_long_cigar_support) {
            EXPECT_FALSE(b.Impl().HasTag("CG"));
        } else {
            EXPECT_TRUE(b.Impl().HasTag("CG"));
        }

        BamWriter writer{LongCigarTests::LongCigarOut, b.header_};
        writer.Write(b);
    }

    {  // read back in
        const auto b = LongCigarTests::ReadLongCigarRecord(LongCigarTests::LongCigarOut);

        EXPECT_EQ(LongCigarTests::numOps, b.CigarData().size());
        if (LongCigarTests::has_native_long_cigar_support) {
            EXPECT_FALSE(b.Impl().HasTag("CG"));
        } else {
            EXPECT_TRUE(b.Impl().HasTag("CG"));
        }
    }
}
