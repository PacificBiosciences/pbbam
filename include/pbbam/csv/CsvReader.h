#ifndef PBBAM_CSV_CSVREADER_H
#define PBBAM_CSV_CSVREADER_H

#include <pbbam/Config.h>

#include <pbbam/TextFileReader.h>
#include <pbbam/csv/CsvTypes.h>
#include <pbbam/internal/QueryBase.h>

#include <filesystem>
#include <string>
#include <vector>

namespace PacBio {
namespace CSV {

///
/// \brief Read CSV/TSV files.
///
/// Supports range-for iteration, returning each CsvRecord (std::map of column name
/// to row value). This avoids the need for clients to have explicit knowledge of
/// column order in the CSV.
///
/// \code
///    CsvReader reader{filename};
///    for (const CsvRecord& record : reader) {
///        const std::string& s = record.at("name");
///        const int x = std::stoi(record.at("age"));
///    }
/// \endcode
///
/// The constructor with no explicit delimiter does an auto-detection step upfront.
/// It will select the best guess among the following {'\t', ',', ' ', '|'}. If
/// the exact delimiter is known for a file, that character may be specified.
/// (This method is preferred, in case the detection heuristic fails.)
///
/// Both plain-text and gzipped input are supported.
///
class CsvReader : public BAM::internal::QueryBase<CsvRecord>
{
public:
    // Auto-detect delimiter
    explicit CsvReader(const std::filesystem::path& filename);

    // Use explicit delimiter. Prefer this ctor if known.
    CsvReader(const std::filesystem::path& filename, char delimiter);

    // Returns list of column names, as seen in file
    const CsvHeader& Header() const;

    // Likely not necessary in client code, range-for is supported.
    bool GetNext(CsvRecord& record) final;

private:
    BAM::TextFileReader reader_;
    CsvHeader header_;
    char delimiter_;
    std::string line_;
    std::vector<std::string> buffer_;
    int lineNumber_ = 0;
};

}  // namespace CSV
}  // namespace PacBio

#endif  // PBBAM_CSV_CSVREADER_H
