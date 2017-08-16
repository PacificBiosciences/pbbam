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

//#define private public

//#include "PbbamTestData.h"
//#include <gtest/gtest.h>
//#include <pbbam/EntireFileQuery.h>

//#include <pbbam/UnmappedReadsQuery.h>
//#include <string>
//using namespace PacBio;
//using namespace PacBio::BAM;
//using namespace std;

//const string inputBamFn1 = PbbamTestsConfig::Data_Dir + "/unmap1.bam";
//const string inputBamFn2 = PbbamTestsConfig::Data_Dir + "/unmap2.bam";

//TEST(UnmappedReadsQueryTest, UnmappedOnlyFile)
//{
//    // open input BAM file
//    BamFile bamFile(inputBamFn1);
//    EXPECT_TRUE(bamFile);

//    // check all records, and save unmapped count
//    int count = 0;
//    int unmappedExpected = 0;
//    EntireFileQuery entireFile(bamFile);
//    EXPECT_TRUE(entireFile);
//    for ( const BamRecord& record : entireFile ) {
//        ++count;
//        if (!record.IsMapped())
//            ++unmappedExpected;
//    }
//    EXPECT_EQ(10, count);
//    EXPECT_EQ(10, unmappedExpected);

//    // query unmapped records only
//    int unmappedObserved = 0;
//    UnmappedReadsQuery unmappedReads(bamFile);
//    EXPECT_TRUE(unmappedReads);
//    for ( const BamRecord& record : unmappedReads ) {
//        EXPECT_FALSE(record.IsMapped());
//        ++unmappedObserved;
//    }
//    EXPECT_EQ(unmappedExpected, unmappedObserved);
//}

//TEST(UnmappedReadsQueryTest, MixedFile)
//{
//    // open input BAM file
//    BamFile bamFile(inputBamFn2);
//    EXPECT_TRUE(bamFile);

//    // check all records, and save unmapped count
//    int count = 0;
//    int unmappedExpected = 0;
//    EntireFileQuery entireFile(bamFile);
//    EXPECT_TRUE(entireFile);
//    for ( const BamRecord& record : entireFile ) {
//        ++count;
//        if (!record.IsMapped())
//            ++unmappedExpected;
//    }
//    EXPECT_EQ(19, count);
//    EXPECT_EQ(9, unmappedExpected);

//    // query unmapped records only
//    int unmappedObserved = 0;
//    UnmappedReadsQuery unmappedReads(bamFile);
//    EXPECT_TRUE(unmappedReads);
//    for ( const BamRecord& record : unmappedReads ) {
//        EXPECT_FALSE(record.IsMapped());
//        ++unmappedObserved;
//    }
//    EXPECT_EQ(unmappedExpected, unmappedObserved);
//}

// TODO: handle no index case

// TODO: additional special cases as needed
