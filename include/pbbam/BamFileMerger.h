// Author: Derek Barnett

#ifndef BAMFILEMERGER_H
#define BAMFILEMERGER_H

#include <pbbam/BamHeader.h>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class DataSet;

class BamFileMerger
{
public:
    /// \brief Runs merger on BAM files.
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] bamFilenames         input filenames
    /// \param[in] outputFilename       resulting BAM output
    /// \param[in] createPbi            if true, creates a PBI alongside output BAM
    /// \param[in] initialOutputHeader  "seed" header (used by pbmerge to set @PG entry)
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const std::vector<std::string>& bamFilenames,
                      const std::string& outputFilename, bool createPbi = true,
                      BamHeader initialOutputHeader = BamHeader{});

    /// \brief Runs merger on a dataset, applying any supplied filters.
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] dataset              provides input filenames & filters
    /// \param[in] outputFilename       resulting BAM output
    /// \param[in] createPbi            if true, creates a PBI alongside output BAM
    /// \param[in] initialOutputHeader  "seed" header (used by pbmerge to set @PG entry)
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const PacBio::BAM::DataSet& dataset, const std::string& outputFilename,
                      bool createPbi = true, BamHeader initialOutputHeader = BamHeader{});
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BAMFILEMERGER_H
