#include "PbbamInternalConfig.h"

#include <pbbam/csv/CsvReader.h>

#include <pbbam/csv/CsvTypes.h>
#include "../ErrnoReason.h"

#include <pbcopper/utility/Ssize.h>

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <array>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PacBio {
namespace CSV {
namespace {

CsvHeader MakeHeader(const std::string& line, const char delimiter)
{
    CsvHeader header;
    boost::split(
        header, line, [delimiter](const char c) { return c == delimiter; },
        boost::token_compress_off);
    return header;
}

char DetectDelimiter(const std::filesystem::path& filename)
{
    constexpr std::array<char, 4> CANDIDATES{'\t', ',', ' ', '|'};
    constexpr int PEEK_LINE_COUNT = 3;

    // grab first few lines
    std::vector<std::string> lines;
    lines.reserve(PEEK_LINE_COUNT);
    {
        std::string line;
        BAM::TextFileReader reader{filename.string()};
        while ((lines.size() < PEEK_LINE_COUNT) && reader.GetNext(line)) {
            lines.push_back(std::move(line));
        }
    }

    //
    // Find delimiters whose number of fields are consistent across each line,
    // including the header. Matching the header should smoke out anything with
    // multi-value columns, like:
    //
    // column_1 column_2 column_3
    //   a        b       1,2,3
    //   c        d       4,5,6
    //
    // Then store the number of fields for any consistent delimiter.
    //
    std::map<char, int> delimiterCounts;
    std::vector<std::string> rowFields;
    std::vector<int> rowFieldCounts;
    for (const char candidate : CANDIDATES) {
        delimiterCounts[candidate] = 0;
        rowFields.clear();
        rowFieldCounts.clear();
        for (const std::string_view line : lines) {
            boost::split(
                rowFields, line, [candidate](const char c) { return c == candidate; },
                boost::token_compress_off);
            rowFieldCounts.push_back(rowFields.size());
        }

        delimiterCounts[candidate] =
            (std::all_of(std::cbegin(rowFieldCounts), std::cend(rowFieldCounts),
                         [&](const int x) { return x == rowFieldCounts.front(); })
                 ? rowFieldCounts.front()
                 : 0);
    }

    // In the (unlikely) case where more than 1 consistent delimter, return the
    // one with the highest number of fields.
    const auto bestCandidateIter =
        std::max_element(std::begin(delimiterCounts), std::end(delimiterCounts),
                         [](const auto& p1, const auto& p2) { return p1.second < p2.second; });
    return bestCandidateIter->first;
}

}  // namespace

CsvReader::CsvReader(const std::filesystem::path& filename)
    : CsvReader{filename, DetectDelimiter(filename)}
{}

CsvReader::CsvReader(const std::filesystem::path& filename, const char delimiter)
    : reader_{filename.string()}, delimiter_{delimiter}
{
    // read 1st line as header
    if (!reader_.GetNext(line_)) {
        std::ostringstream s;
        s << "[pbbam] CSV reader ERROR: could read header:\n"
          << "  file: " << filename;
        BAM::MaybePrintErrnoReason(s);
        throw std::runtime_error{s.str()};
    }
    header_ = MakeHeader(line_, delimiter_);
}

bool CsvReader::GetNext(CsvRecord& record)
{
    // likely EOF
    if (!reader_.GetNext(line_)) {
        return false;
    }
    ++lineNumber_;

    // split & validate
    boost::split(
        buffer_, line_, [this](const char c) { return c == delimiter_; },
        boost::token_compress_off);
    if (buffer_.size() != header_.size()) {
        std::ostringstream s;
        s << "[pbbam] CSV reader ERROR: could not parse record, mismatched column/fields:\n"
          << "      file : " << reader_.Filename() << '\n'
          << "    record : " << lineNumber_ << '\n'
          << "  expected : " << header_.size() << " columns\n"
          << "  observed : " << buffer_.size() << " columns\n";
        throw std::runtime_error{s.str()};
    }

    // make record from line
    for (int i = 0; i < Utility::Ssize(header_); ++i) {
        const std::string& key = header_[i];
        record[key] = std::move(buffer_[i]);
    }

    return true;
}

const CsvHeader& CsvReader::Header() const { return header_; }

}  // namespace CSV
}  // namespace PacBio
