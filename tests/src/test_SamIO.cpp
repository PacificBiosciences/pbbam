#include <pbbam/SamReader.h>
#include <pbbam/SamWriter.h>

#include <cstdint>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/StringUtilities.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_SamReader, can_read_basic_sam)
{
    const std::string bamFilename{PbbamTestsConfig::Data_Dir + "/aligned.bam"};
    const std::string samFilename{PbbamTestsConfig::Data_Dir + "/aligned.sam"};

    std::vector<std::string> bamRecordNames;
    std::vector<std::string> samRecordNames;

    BamReader bamInput{bamFilename};
    for (const auto& b : bamInput) {
        bamRecordNames.push_back(b.FullName());
    }

    SamReader samInput{samFilename};
    for (const auto& b : samInput) {
        std::cout << b.FullName() << '\n';
        samRecordNames.push_back(b.FullName());
    }

    EXPECT_EQ(bamRecordNames, samRecordNames);
}

TEST(BAM_SamReader, handles_zero_byte_file)
{
    try {
        SamReader reader{PbbamTestsConfig::Data_Dir + "/zero_bytes.sam"};
        ASSERT_FALSE("should not get here");
    } catch (const std::exception& e) {
        const std::string msg{e.what()};
        EXPECT_TRUE(msg.find("[pbbam] SAM reader ERROR: could not read from empty input:") !=
                    std::string::npos);
    }
}

TEST(BAM_SamWriter, can_roundtrip_header)
{
    // setup header
    const std::string hdrText{
        "@HD\tVN:1.1\tSO:unknown\tpb:3.0.3\n"
        "@RG\tID:6002b307\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
        "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
        "PU:test\tPM:SEQUEL\n"};

    // write header to file
    const std::string generatedFn = PbbamTestsConfig::GeneratedData_Dir + "/samwriter_hdr_only.sam";
    {
        const BamHeader inputHeader{hdrText};
        SamWriter writer(generatedFn, inputHeader);
    };

    // check header
    {
        std::ifstream f{generatedFn};
        const std::string text{(std::istreambuf_iterator<char>(f)),
                               std::istreambuf_iterator<char>()};
        EXPECT_EQ(hdrText, text);
    }

    // clean up
    std::remove(generatedFn.c_str());
}

TEST(BAM_SamWriter, can_roundtrip_single_record)
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
    record.Impl().CigarData(Data::Cigar{});
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
    tags["zm"] = std::int32_t{100};
    tags["qs"] = std::int32_t{0};
    tags["qe"] = std::int32_t{5};
    tags["np"] = std::int32_t{1};
    tags["rq"] = static_cast<float>(0.6);
    tags["RG"] = std::string{"6002b307"};
    tags["sn"] = std::vector<float>{0.2f, 0.2f, 0.2f, 0.2f};
    record.Impl().Tags(tags);

    const std::string expectedSamRecord{
        "test/100/0_5\t4\t*\t0\t0\t*\t*\t0\t0\tACGTC\t@@@@@\tRG:Z:6002b307\t"
        "np:i:1\tqe:i:5\tqs:i:0\trq:f:0.6\tsn:B:f,0.2,0.2,0.2,0.2\tzm:i:100"};

    // write data to file
    const std::string generatedFn =
        PbbamTestsConfig::GeneratedData_Dir + "/samwriter_hdr_and_record.sam";
    {
        SamWriter writer(generatedFn, inputHeader);
        writer.Write(record);
    };

    // check header & record
    {
        std::ifstream f{generatedFn};
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
    std::remove(generatedFn.c_str());
}

TEST(BAM_SamWriter, can_roundtrip_long_cigar)
{
    const std::string longCigarFn = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/long-cigar-1.7.bam";
    const std::string samFn =
        PacBio::BAM::PbbamTestsConfig::GeneratedData_Dir + "/long-cigar-1.7.sam";

    std::string originalCigar;

    // Generate SAM from long CIGAR BAM
    {
        const BamFile inFile{longCigarFn};
        SamWriter writer{samFn, inFile.Header()};
        EntireFileQuery query{inFile};
        for (auto record : query) {
            originalCigar = record.CigarData().ToStdString();
            writer.Write(record);
        }
    }

    // Verify expected output
    {
        std::ifstream f{samFn};

        std::string line1;
        std::string line2;
        std::string line3;
        std::string line4;

        std::getline(f, line1);
        std::getline(f, line2);
        std::getline(f, line3);
        std::getline(f, line4);

        EXPECT_EQ(0, line1.find("@HD"));
        EXPECT_EQ(0, line2.find("@SQ"));
        EXPECT_EQ(0, line3.find("@PG"));

        // This is _literal_ value stored in the CIGAR field for this long-CIGAR
        // record. The real CIGAR data is stored in the "CG" tag.
        //
        // That literal value does not belong in SAM (as well as the CG tag).
        // CIGAR data should be placed in the standard SAM field.
        //
        EXPECT_EQ(std::string::npos, line4.find("457350S497223N"));
        EXPECT_EQ(std::string::npos, line4.find("CG:B:I,"));

        const auto fields = PacBio::BAM::Split(line4);
        ASSERT_EQ(11, fields.size());
        EXPECT_EQ(originalCigar, fields.at(5));
    }
}
