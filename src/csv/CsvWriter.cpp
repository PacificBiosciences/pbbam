#include "PbbamInternalConfig.h"

#include <pbbam/csv/CsvWriter.h>

#include <boost/algorithm/string.hpp>

#include <sstream>
#include <stdexcept>

namespace PacBio {
namespace CSV {

CsvWriter::CsvWriter(const std::filesystem::path& filename, CsvHeader header, const char delimiter,
                     const std::vector<std::string>& comments)
    : writer_{filename.string()}, header_{std::move(header)}, delimiter_(1, delimiter)
{
    for (const auto& comment : comments) {
        writer_.Write(comment);
    }

    if (header_.empty()) {
        std::ostringstream s;
        s << "[pbbam] CSV writer ERROR: cannot write empty header line\n"
          << "  file: " << filename;
        throw std::runtime_error{s.str()};
    }
    const std::string headerText{boost::join(header_, delimiter_)};
    writer_.Write(headerText);

    fields_.reserve(header_.size());
}

void CsvWriter::Write(const CsvRecord& record)
{
    fields_.clear();
    for (const auto& column : header_) {
        fields_.push_back(record.at(column));
    }

    const std::string recordText{boost::join(fields_, delimiter_)};
    writer_.Write(recordText);
}

}  // namespace CSV
}  // namespace PacBio
