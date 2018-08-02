// Author: Derek Barnett

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/IndexedBamWriter.h>
#include <pbbam/PbiBuilder.h>
#include <pbbam/PbiRawData.h>

TEST(IndexedBamWriter, WritesValidIndex)
{
    using namespace PacBio::BAM;

    const std::string inBam = PbbamTestsConfig::Data_Dir + "/polymerase/internal.subreads.bam";
    const std::string outBam = PbbamTestsConfig::GeneratedData_Dir + "/ibw.bam";
    const std::string outPbi = PbbamTestsConfig::GeneratedData_Dir + "/ibw.bam.pbi";

    const BamFile file{inBam};
    const auto& header = file.Header();

    {  // copy file & generate index

        BamReader reader{file};
        IndexedBamWriter writer{outBam, header};

        BamRecord b;
        while (reader.GetNext(b))
            writer.Write(b);
    }

    // close scope to finalize BAM/PBI output

    {  // check random access using PBI

        const PbiRawData idx{outPbi};
        const auto& offsets = idx.BasicData().fileOffset_;

        BamReader reader{outBam};
        BamRecord b;
        for (size_t i = 0; i < offsets.size(); ++i) {
            auto canRead = [](BamReader& myReader, BamRecord& record,
                              size_t loopI) -> ::testing::AssertionResult {
                if (myReader.GetNext(record))
                    return ::testing::AssertionSuccess() << "i: " << loopI;
                else
                    return ::testing::AssertionFailure() << "i: " << loopI;
            };

            reader.VirtualSeek(offsets.at(i));
            EXPECT_TRUE(canRead(reader, b, i));
        }
    }

    // temp file cleanup
    ::remove(outBam.c_str());
    ::remove(outPbi.c_str());
}
