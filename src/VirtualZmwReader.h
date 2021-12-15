#ifndef PBBAM_VIRTUALZMWREADER_H
#define PBBAM_VIRTUALZMWREADER_H

#include <pbbam/Config.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/virtual/VirtualZmwBamRecord.h>

#include <memory>

namespace PacBio {
namespace BAM {

class VirtualZmwReader
{
public:
    /// \brief Creates a reader that will operate on a primary %BAM file (e.g.
    ///        subread data) and a scraps file, consuming all reads.
    ///
    /// \param[in] primaryBamFilepath hqregion.bam or subreads.bam file path
    /// \param[in] scrapsBamFilepath  scraps.bam file path
    ///
    VirtualZmwReader(const std::string& primaryBamFilepath, const std::string& scrapsBamFilepath);

    /// \brief Creates a reader that will operate on a primary %BAM file (e.g.
    ///        subread data) and a scraps file, respecting the provided PBI
    ///        filter.
    ///
    /// \note All %BAM files must have a corresponding ".pbi" index file to use
    ///       the filter. You may need to call BamFile::EnsurePacBioIndexExists
    ///       before constructing the reader.
    ///
    /// \param[in] primaryBamFilepath hqregion.bam or subreads.bam file path
    /// \param[in] scrapsBamFilepath  scraps.bam file path
    /// \param[in] filter PBI filter criteria
    ///
    VirtualZmwReader(const std::string& primaryBamFilepath, const std::string& scrapsBamFilepath,
                     const PbiFilter& filter);

    VirtualZmwReader() = delete;
    VirtualZmwReader(const VirtualZmwReader&) = delete;
    VirtualZmwReader(VirtualZmwReader&&) = delete;
    VirtualZmwReader& operator=(const VirtualZmwReader&) = delete;
    VirtualZmwReader& operator=(VirtualZmwReader&&) = delete;
    ~VirtualZmwReader();

    /// \returns the BamHeader associated with this reader's "primary" %BAM file
    BamHeader PrimaryHeader() const;

    /// \returns the BamHeader associated with this reader's "scraps" %BAM file
    BamHeader ScrapsHeader() const;

    /// \return the BamHeader associated with the newly stitched BAM data
    BamHeader StitchedHeader() const;

    /// \returns true if more ZMWs are available for reading.
    bool HasNext();

    /// \returns the next stitched polymerase read
    VirtualZmwBamRecord Next();

    /// \returns the next set of reads that belong to one ZMW.
    ///          This enables stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw();

private:
    std::unique_ptr<BamFile> primaryBamFile_;
    std::unique_ptr<BamFile> scrapsBamFile_;
    std::unique_ptr<internal::IQuery> primaryQuery_;
    std::unique_ptr<internal::IQuery> scrapsQuery_;
    internal::IQuery::iterator primaryIt_;
    internal::IQuery::iterator scrapsIt_;
    std::unique_ptr<BamHeader> stitchedHeader_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_VIRTUALZMWREADER_H
