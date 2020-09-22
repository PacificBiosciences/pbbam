// Author: Derek Barnett

#include <pbbam/QualityValues.h>

#include <cstddef>
#include <cstdint>

#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_QualityValue, defaults_to_zero)
{
    const QualityValue value;
    EXPECT_EQ(0, value);
    EXPECT_EQ('!', value.Fastq());
}

TEST(BAM_QualityValue, can_create_from_integer)
{
    const QualityValue zero = 0;
    const QualityValue thirtyThree = 33;
    const QualityValue valid = 42;
    const QualityValue max = 93;
    const QualityValue tooHigh = 94;
    const QualityValue wayTooHigh = std::numeric_limits<int8_t>::max();

    EXPECT_EQ(0, zero);
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

TEST(QualityValueTest, can_create_from_fastq_character)
{
    EXPECT_EQ(0, QualityValue::FromFastq('!'));
    EXPECT_EQ(33, QualityValue::FromFastq('B'));
    EXPECT_EQ(42, QualityValue::FromFastq('K'));
    EXPECT_EQ(93, QualityValue::FromFastq('~'));
}

TEST(BAM_QualityValues, default_is_empty)
{
    const QualityValues qvs;
    EXPECT_TRUE(qvs.empty());
    EXPECT_EQ("", qvs.Fastq());
}

TEST(BAM_QualityValues, can_create_from_integer_vector)
{
    const std::string fastqString{"~~~KKBB!!"};
    const std::vector<uint8_t> values{93, 93, 93, 42, 42, 33, 33, 0, 0};

    QualityValues qvs;
    for (auto qv : values)
        qvs.push_back(qv);
    EXPECT_EQ(fastqString, qvs.Fastq());
}

TEST(BAM_QualityValues, can_create_from_fastq_string)
{
    const std::string fastqString{"~~~KKBB!!"};
    const std::vector<uint8_t> values{93, 93, 93, 42, 42, 33, 33, 0, 0};

    const auto qvs = QualityValues::FromFastq(fastqString);
    EXPECT_EQ(fastqString.size(), qvs.size());
    EXPECT_EQ(values.size(), qvs.size());
    for (size_t i = 0; i < fastqString.size(); ++i) {
        EXPECT_EQ(values.at(i), qvs.at(i));
    }
}
