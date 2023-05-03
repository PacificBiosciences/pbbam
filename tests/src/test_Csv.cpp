#include <pbbam/csv/CsvReader.h>
#include <pbbam/csv/CsvWriter.h>

#include "PbbamTestData.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <tuple>

using namespace PacBio;

namespace CsvReaderTests {

const std::vector<std::string> Comments{
    "#Some comment",
    "##Another comment",
};

// clang-format off
void CheckBasicCsv(CSV::CsvReader& reader) 
{
    // check header
    const CSV::CsvHeader expectedHeader{
        "fruit", "direction", "triforce"
    };
    EXPECT_EQ(expectedHeader, reader.Header());

    // check records
    const std::vector<CSV::CsvRecord> expectedRecords {
        {{"fruit", "apple"},  {"direction", "up"},   {"triforce", "power"}},
        {{"fruit", "banana"}, {"direction", "down"}, {"triforce", "wisdom"}},
        {{"fruit", "orange"}, {"direction", "left"}, {"triforce", "courage"}},
    };

    std::vector<CSV::CsvRecord> observedRecords;
    for (CSV::CsvRecord& record : reader) {
        observedRecords.push_back(std::move(record));
    }
    EXPECT_EQ(observedRecords, expectedRecords);

    // check comments
    EXPECT_EQ(Comments, reader.Comments());
}
// clang-format on

}  // namespace CsvReaderTests

TEST(CSV_CsvReader, throws_on_empty_file)
{
    EXPECT_THROW(CSV::CsvReader reader{BAM::PbbamTestsConfig::Data_Dir + "/csv/empty.csv"},
                 std::runtime_error);
}

TEST(CSV_CsvReader, throws_on_comments_only)
{
    EXPECT_THROW(CSV::CsvReader{BAM::PbbamTestsConfig::Data_Dir + "/csv/comments_only.csv"},
                 std::runtime_error);
}

TEST(CSV_CsvReader, header_only_yields_no_records)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/header_only.csv"};

    CSV::CsvReader reader{fn};
    int count = 0;
    for (const auto& record : reader) {
        std::ignore = record;
        ++count;
    }
    EXPECT_EQ(0, count);

    const CSV::CsvHeader expectedHeader{"fruit", "direction", "triforce"};
    EXPECT_EQ(reader.Header(), expectedHeader);
}

TEST(CSV_CsvReader, can_read_basic_comma_separated)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/comma_separated.csv"};
    CSV::CsvReader reader{fn, ','};
    CsvReaderTests::CheckBasicCsv(reader);
}

TEST(CSV_CsvReader, can_read_gzipped_comma_separated)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/comma_separated.csv.gz"};
    CSV::CsvReader reader{fn, ','};
    CsvReaderTests::CheckBasicCsv(reader);
}

TEST(CSV_CsvReader, can_read_basic_tab_separated)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/tab_separated.csv"};
    CSV::CsvReader reader{fn, '\t'};
    CsvReaderTests::CheckBasicCsv(reader);
}

TEST(CSV_CsvReader, can_read_gzipped_tab_separated)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/tab_separated.csv.gz"};
    CSV::CsvReader reader{fn, '\t'};
    CsvReaderTests::CheckBasicCsv(reader);
}

TEST(CSV_CsvReader, can_handle_internal_comment)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/internal_comment.csv"};
    CSV::CsvReader reader{fn};
    CsvReaderTests::CheckBasicCsv(reader);
}

TEST(CSV_CsvReader, can_handle_end_of_file_comment)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/end_comment.csv"};
    CSV::CsvReader reader{fn};
    CsvReaderTests::CheckBasicCsv(reader);
}

TEST(CSV_CsvReader, throws_if_missing_fields)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/missing_fields.csv"};

    //   | fruit,direction,triforce
    // 1 | apple,up,power
    // 2 | banana,down
    // 3 | orange,left,courage

    try {
        CSV::CsvReader reader{fn, ','};
        for (const auto& record : reader) {
            std::ignore = record;
        }
        ASSERT_FALSE("reached here");
    } catch (const std::exception& e) {
        const std::string_view msg{e.what()};
        EXPECT_TRUE(msg.find("record : 2") != std::string::npos);
        EXPECT_TRUE(msg.find("expected : 3 columns") != std::string::npos);
        EXPECT_TRUE(msg.find("observed : 2 columns") != std::string::npos);
    }
}

TEST(CSV_CsvReader, throws_if_too_many_fields)
{
    const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/too_many_fields.csv"};

    //   | fruit,direction,triforce
    // 1 | apple,up,power
    // 2 | banana,down,wisdom
    // 3 | orange,left,courage,nope

    try {
        CSV::CsvReader reader{fn, ','};
        for (const auto& record : reader) {
            std::ignore = record;
        }
        ASSERT_FALSE("reached here");
    } catch (const std::exception& e) {
        const std::string_view msg{e.what()};
        EXPECT_TRUE(msg.find("record : 3") != std::string::npos);
        EXPECT_TRUE(msg.find("expected : 3 columns") != std::string::npos);
        EXPECT_TRUE(msg.find("observed : 4 columns") != std::string::npos);
    }
}

TEST(CSV_CsvReader, can_autodetect_delimiter)
{
    // tab
    {
        const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir + "/csv/tab_separated.csv"};
        CSV::CsvReader reader{fn};
        CsvReaderTests::CheckBasicCsv(reader);
    }
    // comma
    {
        const std::filesystem::path fn{BAM::PbbamTestsConfig::Data_Dir +
                                       "/csv/comma_separated.csv"};
        CSV::CsvReader reader{fn};
        CsvReaderTests::CheckBasicCsv(reader);
    }
}

TEST(CSV_CsvWriter, can_roundtrip_comma_separated)
{
    const std::filesystem::path inputFn{BAM::PbbamTestsConfig::Data_Dir +
                                        "/csv/comma_separated.csv"};
    const std::filesystem::path outputFn{BAM::PbbamTestsConfig::GeneratedData_Dir +
                                         "/comma_separated_write.csv"};

    {
        CSV::CsvReader reader{inputFn};
        CSV::CsvWriter writer{outputFn, reader.Header(), ',', CsvReaderTests::Comments};
        for (const auto& record : reader) {
            writer.Write(record);
        }
    }
    {
        CSV::CsvReader reader{outputFn};
        CsvReaderTests::CheckBasicCsv(reader);
    }
}

TEST(CSV_CsvWriter, can_roundtrip_gzipped_comma_separated)
{
    const std::filesystem::path inputFn{BAM::PbbamTestsConfig::Data_Dir +
                                        "/csv/comma_separated.csv.gz"};
    const std::filesystem::path outputFn{BAM::PbbamTestsConfig::GeneratedData_Dir +
                                         "/comma_separated_write.csv.gz"};

    {
        CSV::CsvReader reader{inputFn};
        CSV::CsvWriter writer{outputFn, reader.Header(), ',', CsvReaderTests::Comments};
        for (const auto& record : reader) {
            writer.Write(record);
        }
    }
    {
        CSV::CsvReader reader{outputFn};
        CsvReaderTests::CheckBasicCsv(reader);
    }
}

TEST(CSV_CsvWriter, can_roundtrip_tab_separated)
{
    const std::filesystem::path inputFn{BAM::PbbamTestsConfig::Data_Dir + "/csv/tab_separated.csv"};
    const std::filesystem::path outputFn{BAM::PbbamTestsConfig::GeneratedData_Dir +
                                         "/tab_separated_write.csv"};

    {
        CSV::CsvReader reader{inputFn};
        CSV::CsvWriter writer{outputFn, reader.Header(), '\t', CsvReaderTests::Comments};
        for (const auto& record : reader) {
            writer.Write(record);
        }
    }
    {
        CSV::CsvReader reader{outputFn};
        CsvReaderTests::CheckBasicCsv(reader);
    }
}

TEST(CSV_CsvWriter, can_roundtrip_gzipped_tab_separated)
{
    const std::filesystem::path inputFn{BAM::PbbamTestsConfig::Data_Dir +
                                        "/csv/tab_separated.csv.gz"};
    const std::filesystem::path outputFn{BAM::PbbamTestsConfig::GeneratedData_Dir +
                                         "/tab_separated_write.csv.gz"};

    {
        CSV::CsvReader reader{inputFn};
        CSV::CsvWriter writer{outputFn, reader.Header(), '\t', CsvReaderTests::Comments};
        for (const auto& record : reader) {
            writer.Write(record);
        }
    }
    {
        CSV::CsvReader reader{outputFn};
        CsvReaderTests::CheckBasicCsv(reader);
    }
}
