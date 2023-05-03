#ifndef PBBAM_CSV_CSVWRITER_H
#define PBBAM_CSV_CSVWRITER_H

#include <pbbam/Config.h>

#include <pbbam/TextFileWriter.h>
#include <pbbam/csv/CsvTypes.h>

#include <filesystem>
#include <string>
#include <vector>

namespace PacBio {
namespace CSV {

///
/// \brief Write CSV/TSV files.
///
/// Writes CSV records via a std::map<std::string, std::string> (column name, row value).
/// This avoids the need for clients to track column order in the CSV while writing.
/// The ordering is specified by providing the header names in the constructor.
///
/// \code
///    std::vector<std::string> columns{"name", "age"};
///    CsvWriter writer{filename, columns, '\t'};
///    CsvRecord row;
///    for (const auto& player : team) {
///        row["name"] = player.Name;
///        row["age"] = std::to_string(player.Age);
///        writer.Write(row);
//     }
/// \endcode
///
/// Leading comments before the header may be provided. These are written verbatim,
/// so the caller must be sure to provide the desired prefix (typically '#').
///
/// Both plain-text and gzipped output are supported. Use the ".gz" filename suffix
/// to enable compression.
///
class CsvWriter
{
public:
    CsvWriter(const std::filesystem::path& filename, CsvHeader header, char delimiter,
              const std::vector<std::string>& comments = {});
    void Write(const CsvRecord& record);

private:
    BAM::TextFileWriter writer_;
    CsvHeader header_;
    std::string delimiter_;
    std::vector<std::string> fields_;
};

}  // namespace CSV
}  // namespace PacBio

#endif  // PBBAM_CSV_CSVWRITER_H
