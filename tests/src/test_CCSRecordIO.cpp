// Author: Derek Barnett

#include <pbbam/ccs/CCSRecord.h>

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/ccs/CCSRecordFormat.h>
#include <pbbam/ccs/CCSRecordReader.h>
#include <pbbam/ccs/CCSRecordWriter.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::CCS;

// clang-format off

namespace CCSRecordIOTests {

const std::vector<std::string> ValidHeaderText{
    "movie_name=m54238_180925_225123",
    "binding_kit=101-789-500",
    "sequencing_kit=101-789-300",
    "basecaller_version=5.0",
    "framerate=100"
};

const CCSHeader& ValidHeader() {
    static const CCSHeader header{
        "m54238_180925_225123",
        "101-789-500",
        "101-789-300",
        "5.0",
        "100"
    };
    return header;
}

const std::string ValidRecordText{
    "4391137\t0\t459\t2\t0.8\t7.6,13.9,7,12.2\tGATTACA\t13,8,3,14,18,3"
};

const CCSRecord& ValidRecord()
{
    static const CCSRecord record = [](){
        CCSRecord r;
        r.HoleNumber = 4391137;
        r.QueryStart = 0;
        r.QueryEnd = 459;
        r.LocalContextFlags = Data::LocalContextFlags::ADAPTER_AFTER;
        r.Accuracy = 0.8f;
        r.SignalToNoise = {7.6, 13.9, 7.0, 12.2};
        r.Sequence = "GATTACA";
        r.PulseWidths = Data::Frames{std::vector<uint16_t>{13, 8, 3, 14, 18, 3}};
        return r;
    }();
    return record;
}

void CheckHeader(const CCSHeader& expected, const CCSHeader& observed)
{
    EXPECT_EQ(expected.MovieName, observed.MovieName);
    EXPECT_EQ(expected.BindingKit, observed.BindingKit);
    EXPECT_EQ(expected.SequencingKit, observed.SequencingKit);
    EXPECT_EQ(expected.BasecallerVersion, observed.BasecallerVersion);
    EXPECT_EQ(expected.FrameRate, observed.FrameRate);
}

void CheckRecord(const CCSRecord& expected, const CCSRecord& observed)
{
    EXPECT_EQ(expected.HoleNumber, observed.HoleNumber);
    EXPECT_EQ(expected.QueryStart, observed.QueryStart);
    EXPECT_EQ(expected.QueryEnd, observed.QueryEnd);
    EXPECT_EQ(expected.LocalContextFlags, observed.LocalContextFlags);
    EXPECT_EQ(expected.Accuracy, observed.Accuracy);
    EXPECT_EQ(expected.SignalToNoise, observed.SignalToNoise);
    EXPECT_EQ(expected.Sequence, observed.Sequence);

    ASSERT_FALSE(expected.PulseWidths.empty());
    ASSERT_FALSE(observed.PulseWidths.empty());
    EXPECT_TRUE(std::equal(expected.PulseWidths.cbegin(), expected.PulseWidths.cend(),
                           observed.PulseWidths.cbegin()));
}

}  // namespace CCSRecordIOTests

// clang-format on

TEST(CCSRecordIOTest, can_deserialize_valid_header_text)
{
    const auto& lines = CCSRecordIOTests::ValidHeaderText;
    const CCSHeader result = CCSRecordFormat::DeserializeHeader(lines);
    CCSRecordIOTests::CheckHeader(CCSRecordIOTests::ValidHeader(), result);
}

TEST(CCSRecordIOTest, deserialization_throws_on_invalid_header_text)
{
    // clang-format off

    const std::vector<std::string> InvalidHeaderText_Empty{};

    const std::vector<std::string> InvalidHeaderText_EmptyLine{
        "movie_name=m54238_180925_225123=error",
        "",
        "binding_kit=101-789-500",
        "sequencing_kit=101-789-300",
        "basecaller_version=5.0",
        "framerate=100"
    };

    const std::vector<std::string> InvalidHeaderText_ExtraEquals{
        "movie_name=m54238_180925_225123=error",
        "binding_kit=101-789-500",
        "sequencing_kit=101-789-300",
        "basecaller_version=5.0",
        "framerate=100"
    };

    const std::vector<std::string> InvalidHeaderText_MissingEquals{
        "movie_name=m54238_180925_225123",
        "binding_kit101-789-500",
        "sequencing_kit=101-789-300",
        "basecaller_version=5.0",
        "framerate=100"
    };

    const std::vector<std::string> InvalidHeaderText_UnknownFieldName{
        "movie_name=m54238_180925_225123",
        "binding_kit101-789-500",
        "sequencing_kit=101-789-300",
        "basecaller_version=5.0",
        "framerate=100",
        "this=does_not_exist"
    };


    EXPECT_THROW(CCSRecordFormat::DeserializeHeader(InvalidHeaderText_Empty),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeHeader(InvalidHeaderText_EmptyLine),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeHeader(InvalidHeaderText_ExtraEquals),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeHeader(InvalidHeaderText_MissingEquals),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeHeader(InvalidHeaderText_UnknownFieldName),
                 std::runtime_error);

    // clang-format on
}

TEST(CCSRecordIOTest, can_serialize_header)
{
    const auto& expected = CCSRecordIOTests::ValidHeaderText;
    const auto lines = CCSRecordFormat::SerializeHeader(CCSRecordIOTests::ValidHeader());
    EXPECT_TRUE(std::equal(lines.cbegin(), lines.cend(), expected.cbegin()));
}

TEST(CCSRecordIOTest, can_deserialize_valid_record)
{
    const auto& line = CCSRecordIOTests::ValidRecordText;
    const auto observed = CCSRecordFormat::DeserializeRecord(line);
    CCSRecordIOTests::CheckRecord(CCSRecordIOTests::ValidRecord(), observed);
}

TEST(CCSRecordIOTest, deserialization_throws_on_invalid_record)
{
    // clang-format off

    const std::string InvalidRecordText_Empty;

    const std::string InvalidRecordText_TooFewFields{"4391137\t0\t459\t2"};

    const std::string InvalidRecordText_TooManyFields{
        "4391137\t0\t459\t2\t0.8\t7.6,13.9,7,12.2\tGATTACA\t13,8,3,14,18,3\ttoo\tmany\fields"};

    const std::string InvalidRecordText_WrongFieldDelmiter{
        "4391137 0 459 2 0.8 7.6,13.9,7,12.2 GATTACA 13,8,3,14,18,3"};

    const std::string InvalidRecordText_WrongSnrDelmiter{
        "4391137\t0\t459\t2\t0.8\t7.6-13.9-7-12.2\tGATTACA\t13,8,3,14,18,3"};

    EXPECT_THROW(CCSRecordFormat::DeserializeRecord(InvalidRecordText_Empty),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeRecord(InvalidRecordText_TooFewFields),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeRecord(InvalidRecordText_TooManyFields),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeRecord(InvalidRecordText_WrongFieldDelmiter),
                 std::runtime_error);
    EXPECT_THROW(CCSRecordFormat::DeserializeRecord(InvalidRecordText_WrongSnrDelmiter),
                 std::runtime_error);

    // clang-format on
}

TEST(CCSRecordIOTest, can_serialize_record)
{
    const auto& expected = CCSRecordIOTests::ValidRecordText;
    const auto result = CCSRecordFormat::SerializeRecord(CCSRecordIOTests::ValidRecord());
    EXPECT_EQ(expected, result);
}

TEST(CCSRecordIOTest, can_do_round_trip_read_and_write_to_iostreams)
{
    const size_t NumOutputRecords = 3;

    // write to ostream
    std::ostringstream output;
    {
        CCSRecordWriter writer{CCSRecordIOTests::ValidHeader(), output};
        for (size_t i = 0; i < NumOutputRecords; ++i) {
            writer.Write(CCSRecordIOTests::ValidRecord());
        }
    }

    // use ostream contents as istream
    std::istringstream input;
    input.str(output.str());

    // check contents
    CCSRecordReader reader{input};
    CCSRecordIOTests::CheckHeader(CCSRecordIOTests::ValidHeader(), reader.Header());

    size_t recordCount = 0;
    for (const auto& record : reader) {
        CCSRecordIOTests::CheckRecord(CCSRecordIOTests::ValidRecord(), record);
        ++recordCount;
    }
    EXPECT_EQ(NumOutputRecords, recordCount);
}

TEST(CCSRecordTest, can_convert_to_read)
{
    const int32_t holeNumber = 77;
    const Data::Position qStart = 1000;
    const Data::Position qEnd = 1010;
    const Data::LocalContextFlags ctxtFlags =
        Data::LocalContextFlags::ADAPTER_BEFORE | Data::LocalContextFlags::ADAPTER_AFTER;
    const Data::Accuracy acc = 0.95f;
    const Data::SNR snr{0.4, 0.4, 0.4, 0.4};
    const std::string seq{"GGTTAACCAA"};
    const Data::Frames pw{3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    const std::string movie = "movie";
    const std::string chemistry = "chemistry";

    const CCS::CCSRecord ccsRecord{holeNumber, qStart, qEnd, ctxtFlags, acc, snr, seq, pw};

    const auto read = ccsRecord.ToRead(movie, chemistry);
    EXPECT_EQ(holeNumber, read.Id.HoleNumber);
    EXPECT_EQ(qStart, read.QueryStart);
    EXPECT_EQ(qEnd, read.QueryEnd);
    EXPECT_EQ(ctxtFlags, read.Flags);
    EXPECT_EQ(acc, read.ReadAccuracy);
    EXPECT_EQ(snr, read.SignalToNoise);
    EXPECT_EQ(seq, read.Seq);
    EXPECT_EQ(pw, read.PulseWidth);
}
