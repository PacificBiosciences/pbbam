#ifndef PBBAM_ZMWREADSTITCHER_H
#define PBBAM_ZMWREADSTITCHER_H

#include <pbbam/Config.h>

#include <memory>
#include <string>
#include <vector>

#include <pbbam/BamRecord.h>
#include <pbbam/DataSet.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/virtual/VirtualZmwBamRecord.h>

namespace PacBio {
namespace BAM {

/// \brief The ZmwReadStitcher class provides an interface for re-stitching
///        "virtual" polymerase reads from their constituent parts.
///
/// \note This reader requires that any input %BAM files also have associated PBI
///       files available for query. See BamFile::EnsurePacBioIndexExists .
///
class PBBAM_EXPORT ZmwReadStitcher
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// entire file, from BAM names
    ZmwReadStitcher(std::string primaryBamFilePath, std::string scrapsBamFilePath);

    /// filtered input from BAM names
    ZmwReadStitcher(std::string primaryBamFilePath, std::string scrapsBamFilePath,
                    PbiFilter filter);

    /// maybe filtered, from DataSet input
    ZmwReadStitcher(const DataSet& dataset);

    ZmwReadStitcher(ZmwReadStitcher&&) noexcept;
    ZmwReadStitcher& operator=(ZmwReadStitcher&&) noexcept;
    ~ZmwReadStitcher();

    /// \}

    /// \name File Headers
    /// \{

    /// \returns the BamHeader associated with this reader's "primary" %BAM file
    BamHeader PrimaryHeader() const;

    /// \returns the BamHeader associated with this reader's "scraps" %BAM file
    BamHeader ScrapsHeader() const;

    /// \return the BamHeader associated with the newly stitched BAM data
    BamHeader StitchedHeader() const;

    /// \}

    /// \name Stitched Record Reading
    ///

    /// \returns true if more ZMWs are available for reading.
    bool HasNext();

    /// \returns the next stitched polymerase read
    VirtualZmwBamRecord Next();

    /// \returns the next set of reads that belong to one ZMW.
    ///          This enables stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw();

    /// \}

private:
    class ZmwReadStitcherPrivate;
    std::unique_ptr<ZmwReadStitcherPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWREADSTITCHER_H
