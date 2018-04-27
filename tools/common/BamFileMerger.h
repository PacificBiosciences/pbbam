// Author: Derek Barnett

#ifndef BAMFILEMERGER_H
#define BAMFILEMERGER_H

#include <pbbam/DataSet.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/ProgramInfo.h>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {
namespace common {

class BamFileMerger
{
public:
    /// \brief Runs merger on a dataset, applying any supplied filters.
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] dataset          provides input filenames & filters
    /// \param[in] outputFilename   resulting BAM output
    /// \param[in] mergeProgram     info about the calling program. Adds a @PG entry to merged header.
    /// \param[in] createPbi        if true, creates a PBI alongside output BAM
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const PacBio::BAM::DataSet& dataset, const std::string& outputFilename,
                      const PacBio::BAM::ProgramInfo& mergeProgram = PacBio::BAM::ProgramInfo(),
                      bool createPbi = true);
};

}  // namespace common
}  // namespace BAM
}  // namespace PacBio

#include "BamFileMerger.inl"

#endif  // BAMFILEMERGER_H
