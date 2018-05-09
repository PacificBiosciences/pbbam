// Author: Derek Barnett

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/SamWriter.h>
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(SamWriterTest, HeaderOk)
{
    // setup header
    const std::string hdrText{
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.3\n"
        "@RG\tID:6002b307\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
        "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
        "PU:test\tPM:SEQUEL\n"};

    EXPECT_NO_THROW({
        // write header to file
        const std::string generatedFn =
            PbbamTestsConfig::GeneratedData_Dir + "/samwriter_hdr_only.sam";
        {
            const BamHeader inputHeader(hdrText);
            SamWriter writer(generatedFn, inputHeader);
            //            ()writer;
        };

        // check header
        {
            std::ifstream f(generatedFn);
            const std::string text((std::istreambuf_iterator<char>(f)),
                                   std::istreambuf_iterator<char>());
            EXPECT_EQ(hdrText, text);
        }

        // clean up
        remove(generatedFn.c_str());
    });
}

TEST(SamWriterTest, SingleRecordOk)
{

    // setup header
    const std::string hdrLine1{"@HD\tVN:1.1\tSO:unknown\tpb:3.0.3"};
    const std::string hdrLine2{
        "@RG\tID:6002b307\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
        "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
        "PU:test\tPM:SEQUEL"};
    const std::string hdrText = hdrLine1 + "\n" + hdrLine2 + "\n";
    const BamHeader inputHeader(hdrText);

    // setup record
    BamRecord record(inputHeader);
    record.Impl().Name("test/100/0_5");
    record.Impl().SetSequenceAndQualities("ACGTC", 5, "@@@@@");
    record.Impl().CigarData("");
    record.Impl().Bin(0);
    record.Impl().Flag(0);
    record.Impl().InsertSize(0);
    record.Impl().MapQuality(0);
    record.Impl().MatePosition(-1);
    record.Impl().MateReferenceId(-1);
    record.Impl().Position(-1);
    record.Impl().ReferenceId(-1);
    record.Impl().SetMapped(false);

    TagCollection tags;
    tags["zm"] = int32_t{100};
    tags["qs"] = int32_t{0};
    tags["qe"] = int32_t{5};
    tags["np"] = int32_t{1};
    tags["rq"] = static_cast<float>(0.6);
    tags["RG"] = std::string{"6002b307"};
    tags["sn"] = std::vector<float>{0.2f, 0.2f, 0.2f, 0.2f};
    record.Impl().Tags(tags);

    const std::string expectedSamRecord{
        "test/100/0_5\t4\t*\t0\t0\t*\t*\t0\t0\tACGTC\t@@@@@\tRG:Z:6002b307\t"
        "np:i:1\tqe:i:5\tqs:i:0\trq:f:0.6\tsn:B:f,0.2,0.2,0.2,0.2\tzm:i:100"};

    EXPECT_NO_THROW({
        // write data to file
        const std::string generatedFn =
            PbbamTestsConfig::GeneratedData_Dir + "/samwriter_hdr_and_record.sam";
        {
            SamWriter writer(generatedFn, inputHeader);
            writer.Write(record);
        };

        // check header & record
        {
            std::ifstream f(generatedFn);
            std::string line1;
            std::string line2;
            std::string line3;
            std::getline(f, line1);
            std::getline(f, line2);
            std::getline(f, line3);
            EXPECT_EQ(hdrLine1, line1);
            EXPECT_EQ(hdrLine2, line2);
            EXPECT_EQ(expectedSamRecord, line3);
        }

        // cleanup
        remove(generatedFn.c_str());
    });
}
