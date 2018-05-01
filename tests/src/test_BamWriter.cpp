// Author: Derek Barnett

#include <gtest/gtest.h>
#include <cstdint>
#include <string>
#include <vector>

#include "PbbamTestData.h"

#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BamWriterTest, SingleWrite_UserRecord)
{
    const std::string fullName = "test/100/0_5";
    const std::string rgId = "6002b307";
    const std::vector<float> expectedSnr = {0.2, 0.2, 0.2, 0.2};

    // setup header
    const std::string hdrText = {
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.1\n"
        "@RG\tID:6002b307\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
        "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
        "PU:test\tPM:SEQUEL\n"};
    BamHeader inputHeader(hdrText);

    // setup record
    BamRecord bamRecord(inputHeader);
    bamRecord.Impl().Name(fullName);
    bamRecord.Impl().SetSequenceAndQualities("ACGTC", 5);
    bamRecord.Impl().CigarData("");
    bamRecord.Impl().Bin(0);
    bamRecord.Impl().Flag(0);
    bamRecord.Impl().InsertSize(0);
    bamRecord.Impl().MapQuality(0);
    bamRecord.Impl().MatePosition(-1);
    bamRecord.Impl().MateReferenceId(-1);
    bamRecord.Impl().Position(-1);
    bamRecord.Impl().ReferenceId(-1);
    bamRecord.Impl().SetMapped(false);

    TagCollection tags;
    tags["zm"] = int32_t{100};
    tags["qs"] = int32_t{0};
    tags["qe"] = int32_t{5};
    tags["np"] = int32_t{1};
    tags["rq"] = static_cast<float>(0.6);
    tags["RG"] = rgId;
    tags["sn"] = expectedSnr;
    bamRecord.Impl().Tags(tags);

    // write record to file
    const std::string generatedBamFn =
        PbbamTestsConfig::GeneratedData_Dir + "/bamwriter_generated.bam";
    {
        BamWriter writer(generatedBamFn, inputHeader);
        writer.Write(bamRecord);
    }

    // check written header
    BamFile file(generatedBamFn);
    const auto header = file.Header();
    EXPECT_EQ(std::string("1.1"), header.Version());
    EXPECT_EQ(std::string("unknown"), header.SortOrder());
    EXPECT_EQ(std::string("3.0.1"), header.PacBioBamVersion());

    // check written record
    EntireFileQuery entireFile(file);
    auto firstIter = entireFile.begin();
    auto record = *firstIter;
    EXPECT_EQ(std::string("ACGTC"), record.Sequence());
    EXPECT_EQ(std::string("test/100/0_5"), record.FullName());
    EXPECT_TRUE(record.HasHoleNumber());
    EXPECT_TRUE(record.HasNumPasses());
    EXPECT_TRUE(record.HasQueryEnd());
    EXPECT_TRUE(record.HasQueryStart());
    EXPECT_TRUE(record.HasReadAccuracy());
    EXPECT_TRUE(record.HasSignalToNoise());
    EXPECT_EQ(100, record.HoleNumber());
    EXPECT_EQ(1, record.NumPasses());
    EXPECT_EQ(0, record.QueryStart());
    EXPECT_EQ(5, record.QueryEnd());
    EXPECT_EQ(expectedSnr, record.SignalToNoise());
    EXPECT_EQ(rgId, record.ReadGroupId());

    // clean up
    remove(generatedBamFn.c_str());
}
