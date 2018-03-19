// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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
            auto canRead = [](BamReader& reader, BamRecord& record,
                              size_t i) -> ::testing::AssertionResult {
                if (reader.GetNext(record))
                    return ::testing::AssertionSuccess() << "i: " << i;
                else
                    return ::testing::AssertionFailure() << "i: " << i;
            };

            reader.VirtualSeek(offsets.at(i));
            EXPECT_TRUE(canRead(reader, b, i));
        }
    }

    // temp file cleanup
    ::remove(outBam.c_str());
    ::remove(outPbi.c_str());
}
