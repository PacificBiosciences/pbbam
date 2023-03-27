#include "PbbamInternalConfig.h"

#include <pbbam/csv/CsvWriter.h>

#include <boost/algorithm/string.hpp>

namespace PacBio {
namespace CSV {

CsvWriter::CsvWriter(const std::filesystem::path& filename, CsvHeader header, const char delimiter)
    : writer_{filename.string()}, header_{std::move(header)}, delimiter_(1, delimiter)
{
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
