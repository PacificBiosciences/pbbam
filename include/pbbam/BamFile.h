// File Description
/// \file BamFile.h
/// \brief Defines the BamFile class.
//
// Author: Derek Barnett

#ifndef BAMFILE_H
#define BAMFILE_H

#include <cstdint>
#include <string>

#include "pbbam/BamHeader.h"
#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

namespace internal {
class BamFilePrivate;
}

/// \brief The BamFile class represents a %BAM file.
///
/// It provides access to header metadata and methods for finding/creating
/// associated index files.
///
class PBBAM_EXPORT BamFile
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a BamFile object on the provided \p filename &
    ///        loads header information.
    ///
    /// \param[in] filename %BAM filename
    /// \throws std::exception on failure to open %BAM file for reading
    ///
    BamFile(std::string filename);

    BamFile(const BamFile& other);
    BamFile(BamFile&& other);
    BamFile& operator=(const BamFile& other);
    BamFile& operator=(BamFile&& other);
    ~BamFile();

    /// \}

public:
    /// \name Index & Filename Methods
    /// \{

    /// \brief Creates a ".pbi" file for this %BAM file.
    ///
    /// \note Existing index file will be overwritten. Use
    ///       EnsurePacBioIndexExists() if this is not desired.
    ///
    /// \throws if PBI file could not be properly created and/or
    ///         written to disk
    ///
    void CreatePacBioIndex() const;

    /// \brief Creates a ".bai" file for this %BAM file.
    ///
    /// \note Existing index file will be overwritten. Use
    ///       EnsureStandardIndexExists() if this is not desired.
    ///
    /// \throws if BAI file could not be properly created (e.g. this
    ///         %BAM is not coordinate-sorted) or could not be written to disk
    ///
    void CreateStandardIndex() const;

    /// \brief Creates a ".pbi" file if one does not exist or is older than its
    ///        %BAM file.
    ///
    /// Equivalent to:
    /// \code{.cpp}
    ///    if (!file.PacBioIndexExists())
    ///        file.CreatePacBioIndex();
    /// \endcode
    ///
    /// \note As of v0.4.02+, no timestamp check is performed. Previously we requr
    /// with an additional timestamp check.
    ///
    /// \throws if PBI file could not be properly created and/or
    ///         written to disk
    ///
    void EnsurePacBioIndexExists() const;

    /// \brief Creates a ".bai" file if one does not exist or is older than its
    ///        %BAM file.
    ///
    /// Equivalent to:
    /// \code{.cpp}
    ///    if (!file.StandardIndexExists())
    ///        file.CreateStandardIndex();
    /// \endcode
    ///
    /// \note As of v0.4.2, no timestamp check is performed.
    ///
    /// \throws if BAI file could not be properly created (e.g. this
    ///         %BAM is not coordinate-sorted) or could not be written to disk
    ///
    void EnsureStandardIndexExists() const;

    /// \returns %BAM filename
    const std::string& Filename() const;

    /// \returns true if %BAM file has EOF marker (empty BGZF block). Streamed
    ///          input (filename: "-")
    bool HasEOF() const;

    /// \returns true if ".pbi" exists and is newer than this %BAM file.
    bool PacBioIndexExists() const;

    /// \returns filename of %PacBio index file (".pbi")
    /// \note No guarantee is made on the existence of this file.
    ///       This method simply returns the expected filename.
    std::string PacBioIndexFilename() const;

    /// \returns true if ".pbi" has a more recent timestamp than this file
    bool PacBioIndexIsNewer() const;

    /// \returns true if ".bai" exists
    bool StandardIndexExists() const;

    /// \note No guarantee is made on the existence of this file.
    ///       This method simply returns the expected filename.
    std::string StandardIndexFilename() const;

    /// \returns true if ".bai" has a more recent timestamp than this file
    bool StandardIndexIsNewer() const;

    /// \}

public:
    /// \name File Header Data
    /// \{

    /// \returns true if header metadata has this reference name
    bool HasReference(const std::string& name) const;

    /// \returns const reference to BamHeader containing the file's metadata
    const BamHeader& Header() const;

    /// \returns true if file is a %PacBio %BAM file (i.e. has non-empty version
    ///          associated with header "pb" tag)
    bool IsPacBioBAM() const;

    /// \returns ID for reference \p name (can be used for e.g.
    ///          GenomicIntervalQuery), or -1 if not found
    int ReferenceId(const std::string& name) const;

    /// \return name of reference matching \p id, empty string if not found
    std::string ReferenceName(const int id) const;

    /// \returns length of requested reference \p name. 0 if not found
    uint32_t ReferenceLength(const std::string& name) const;

    /// \returns length of requested reference \p id. 0 if not found
    uint32_t ReferenceLength(const int id) const;

    /// \}

public:
    /// \name Additional Attributes
    /// \{

    /// \returns virtual offset of first alignment. Intended mostly for internal
    ///          use. Note that this is a BGZF \b virtual offset, not a
    ///          'normal' file position.
    ///
    int64_t FirstAlignmentOffset() const;

    /// \}

private:
    std::unique_ptr<internal::BamFilePrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BAMFILE_H
