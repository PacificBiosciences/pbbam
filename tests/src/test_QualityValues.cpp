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

#include <cstddef>
#include <cstdint>
#include <limits>

#include <gtest/gtest.h>

#define private public

#include <pbbam/QualityValues.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(QualityValueTest, DefaultsOk)
{
    const QualityValue value;
    EXPECT_EQ(0,   value);
    EXPECT_EQ('!', value.Fastq());
}

TEST(QualityValueTest, FromNumber)
{
    const QualityValue zero(0);
    const QualityValue thirtyThree(33);
    const QualityValue valid(42);
    const QualityValue max(93);
    const QualityValue tooHigh(94);
    const QualityValue wayTooHigh(std::numeric_limits<int8_t>::max());

    EXPECT_EQ(0,  zero);
    EXPECT_EQ(33, thirtyThree);
    EXPECT_EQ(42, valid);
    EXPECT_EQ(93, max);
    EXPECT_EQ(93, tooHigh);
    EXPECT_EQ(93, wayTooHigh);

    EXPECT_EQ('!', zero.Fastq());
    EXPECT_EQ('B', thirtyThree.Fastq());
    EXPECT_EQ('K', valid.Fastq());
    EXPECT_EQ('~', max.Fastq());
    EXPECT_EQ('~', tooHigh.Fastq());
    EXPECT_EQ('~', wayTooHigh.Fastq());
}

TEST(QualityValueTest, FromFastq)
{
    const QualityValue zero        = QualityValue::FromFastq('!');
    const QualityValue thirtyThree = QualityValue::FromFastq('B');
    const QualityValue valid       = QualityValue::FromFastq('K');
    const QualityValue max         = QualityValue::FromFastq('~');

    EXPECT_EQ(0,  zero);
    EXPECT_EQ(33, thirtyThree);
    EXPECT_EQ(42, valid);
    EXPECT_EQ(93, max);
}

TEST(QualityValuesTest, Default)
{
    const QualityValues qvs;
    EXPECT_TRUE(qvs.empty());
    EXPECT_EQ(string(), qvs.Fastq());
}

TEST(QualityValuesTest, FromNumbers)
{
    const string fastqString = "~~~KKBB!!";
    const vector<uint8_t> values = { 93, 93, 93, 42, 42, 33, 33, 0, 0 };

    QualityValues qvs;
    for (auto qv : values)
        qvs.push_back(qv);
    EXPECT_EQ(fastqString, qvs.Fastq());
}

TEST(QualityValuesTest, FromFastq)
{
    const string fastqString = "~~~KKBB!!";
    const vector<uint8_t> values = { 93, 93, 93, 42, 42, 33, 33, 0, 0 };

    const QualityValues& qvs = QualityValues::FromFastq(fastqString);
    EXPECT_EQ(fastqString.size(), qvs.size());
    EXPECT_EQ(values.size(), qvs.size());
    for (size_t i = 0; i < fastqString.size(); ++i)
        EXPECT_EQ(values.at(i), qvs.at(i));
}
