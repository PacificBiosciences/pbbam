#ifndef PBBAM_BAMFILEMERGER_H
#define PBBAM_BAMFILEMERGER_H

#include <pbbam/Config.h>

#include <pbbam/DataSet.h>
#include <pbbam/IRecordWriter.h>
#include <pbbam/ProgramInfo.h>

#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

struct BamFileMerger
{
    /// \brief Runs merger on BAM files.
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] bamFilenames     input filenames
    /// \param[in] outputFilename   resulting BAM output
    /// \param[in] createPbi        if true, creates a PBI alongside output BAM
    /// \param[in] pgInfo           allows client applications to add its @PG entry to merged header
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const std::vector<std::string>& bamFilenames,
                      const std::string& outputFilename, bool createPbi = true,
                      const ProgramInfo& pgInfo = ProgramInfo{});

    /// \brief Runs merger on a dataset (applying any filters).
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] dataset          provides input filenames & filters
    /// \param[in] outputFilename   resulting BAM output
    /// \param[in] createPbi        if true, creates a PBI alongside output BAM
    /// \param[in] pgInfo           allows client applications to add its @PG entry to merged header
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const DataSet& dataset, const std::string& outputFilename,
                      bool createPbi = true, const ProgramInfo& pgInfo = ProgramInfo{});

    /// \brief Runs merger on BAM files, writing to provided writer.
    ///
    /// \param[in] bamFilenames     input filenames
    /// \param[in] writer           output destination
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const std::vector<std::string>& bamFilenames, IRecordWriter& writer);

    /// \brief Runs merger on a dataset (applying any filters), writing to
    ///        provided writer.
    ///
    /// \param[in] dataset      provides input filenames & filters
    /// \param[in] writer       output destination
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const DataSet& dataset, IRecordWriter& writer);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BAMFILEMERGER_H
