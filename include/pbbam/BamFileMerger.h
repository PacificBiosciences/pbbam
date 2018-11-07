// Author: Derek Barnett

#ifndef BAMFILEMERGER_H
#define BAMFILEMERGER_H

#include <pbbam/ProgramInfo.h>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class DataSet;
class IRecordWriter;

class BamFileMerger
{
public:
    /// \brief Runs merger on BAM files.
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] bamFilenames      input filenames
    /// \param[in] outputFilename    resulting BAM output
    /// \param[in] createPbi         if true, creates a PBI alongside output BAM
    /// \param[in] pgInfo            allows client applications to add its @PG entry to merged header
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const std::vector<std::string>& bamFilenames,
                      const std::string& outputFilename, bool createPbi = true,
                      const ProgramInfo& pgInfo = ProgramInfo{});

    /// \brief Runs merger on a dataset, applying any supplied filters.
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] dataset              provides input filenames & filters
    /// \param[in] outputFilename       resulting BAM output
    /// \param[in] createPbi            if true, creates a PBI alongside output BAM
    /// \param[in] pgInfo            allows client applications to add its @PG entry to merged header
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const PacBio::BAM::DataSet& dataset, const std::string& outputFilename,
                      bool createPbi = true, const ProgramInfo& pgInfo = ProgramInfo{});

    static void Merge(const std::vector<std::string>& bamFilenames, IRecordWriter& writer);

    static void Merge(const PacBio::BAM::DataSet& dataset, IRecordWriter& writer);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BAMFILEMERGER_H
