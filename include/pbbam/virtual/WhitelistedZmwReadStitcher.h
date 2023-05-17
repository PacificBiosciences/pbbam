#ifndef PBBAM_WHITELISTEDZMWREADSTITCHER_H
#define PBBAM_WHITELISTEDZMWREADSTITCHER_H

#include <pbbam/Config.h>

#include <pbbam/BamRecord.h>
#include <pbbam/DataSet.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/virtual/VirtualZmwBamRecord.h>

#include <memory>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief The WhitelistedZmwReadStitcher class provides an interface for
///        re-stitching "virtual" ZMW reads from their constituent parts,
///        limiting results to only those reads originating from a 'whitelist'
///         of ZMW hole numbers.
///
/// Whitelisted ZMWs that are not present in both primary and scraps BAMs
/// will be "pre-removed." This ensures that, given client code like this:
///
/// \include code/WhitelistedZmwReadStitcher.txt
///
/// each iteration will always provide valid data - either a valid virtual
/// record from Next() or a non-empty vector from NextRaw().
///
/// \note This reader requires that both input %BAM files also have associated
///       PBI files available for query. See BamFile::EnsurePacBioIndexExists .
///
class PBBAM_EXPORT WhitelistedZmwReadStitcher
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a reader that will operate on a primary %BAM file (e.g. subread data)
    ///        and a scraps file, using a ZMW whitelist to filter the input.
    ///
    /// \param[in] zmwWhitelist         list of ZMWs to restrict iteration over
    /// \param[in] primaryBamFilePath   hqregion.bam or subreads.bam file path
    /// \param[in] scrapsBamFilePath    scraps.bam file path
    ///
    /// \note This reader requires that both input %BAM files also have associated PBI
    ///       files available for query. See BamFile::EnsurePacBioIndexExists .
    ///
    /// \throws std::runtime_error if any files (*.bam and/or *.pbi) were not available for reading, or
    ///         if malformed data encountered
    ///
    WhitelistedZmwReadStitcher(const std::vector<std::int32_t>& zmwWhitelist,
                               const std::string& primaryBamFilePath,
                               const std::string& scrapsBamFilePath);

    WhitelistedZmwReadStitcher(WhitelistedZmwReadStitcher&&) noexcept;
    WhitelistedZmwReadStitcher& operator=(WhitelistedZmwReadStitcher&&) noexcept;
    ~WhitelistedZmwReadStitcher();

    /// \}

    /// \name Stitched Record Reading
    /// \{

    /// \returns true if more ZMWs are available for reading.
    bool HasNext() const;

    /// \returns the re-stitched polymerase read from the next ZMW in the whitelist
    VirtualZmwBamRecord Next();

    /// \returns the set of reads that belong to the next ZMW in the whitelist.
    ///          This enables stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw();

    /// \}

    /// \name File Headers
    /// \{

    /// \returns the BamHeader associated with this reader's "primary" %BAM file
    BamHeader PrimaryHeader() const;

    /// \returns the BamHeader associated with this reader's "scraps" %BAM file
    BamHeader ScrapsHeader() const;

    /// \}

private:
    class WhitelistedZmwReadStitcherPrivate;
    std::unique_ptr<WhitelistedZmwReadStitcherPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_WHITELISTEDZMWREADSTITCHER
