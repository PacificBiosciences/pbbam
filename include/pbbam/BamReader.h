// File Description
/// \file BamReader.h
/// \brief Defines the BamReader class.
//
// Author: Derek Barnett

#ifndef BAMREADER_H
#define BAMREADER_H

#include <cstdint>
#include <memory>
#include <string>

#include <htslib/sam.h>

#include "pbbam/BamFile.h"
#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/GenomicInterval.h"

namespace PacBio {
namespace BAM {

namespace internal {
struct BamReaderPrivate;
}

/// \brief The BamReader class provides basic read-access to a %BAM file.
///
/// The base-class implementation provides a sequential read-through of BAM
/// records. Derived classes may implement other access schemes (e.g. genomic
/// region, PBI-enabled record filtering).
///
class PBBAM_EXPORT BamReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Opens BAM file for reading.
    ///
    /// \param[in] fn %BAM filename
    /// \throws std::runtime_error if failed to open
    ///
    explicit BamReader(std::string fn);

    /// \brief Opens BAM file for reading.
    ///
    /// \param[in] bamFile BamFile object
    /// \throws std::runtime_error if failed to open
    ///
    explicit BamReader(BamFile bamFile);

    virtual ~BamReader();

    /// \}

public:
    /// \name BAM File Attributes
    /// \{

    /// \returns the underlying BamFile
    const BamFile& File() const;

    /// \returns %BAM filename
    const std::string& Filename() const;

    /// \returns BamHeader object from %BAM header contents
    const BamHeader& Header() const;

    /// \}

public:
    /// \name BAM File I/O
    /// \{

    /// \brief Fetches the "next" %BAM record.
    ///
    /// Default implementation will read records until EOF. Derived readers may
    /// use additional criteria to decide which record is "next" and when
    /// reading is done.
    ///
    /// \param[out] record  next BamRecord object. Should not be used if method
    ///                     returns false.
    ///
    /// \returns true if record was read successfully. Returns false if EOF (or
    ///          end of iterator in derived readers). False is not an error,
    ///          it indicates "end of data".
    ///
    /// \throws std::runtime_error if failed to read from file (e.g. possible
    ///         truncated or corrupted file).
    ///
    bool GetNext(BamRecord& record);

    /// \brief Seeks to virtual offset in %BAM.
    ///
    /// \note This is \b NOT a normal file offset, but the virtual offset used
    ///       in %BAM indexing.
    ///
    /// \throws std::runtime_error if failed to seek
    ///
    void VirtualSeek(int64_t virtualOffset);

    /// \returns current (virtual) file position.
    ///
    /// \note This is \b NOT a normal file offset, but the virtual offset used
    ///       in %BAM indexing.
    ///
    int64_t VirtualTell() const;

    /// \}

protected:
    /// \name BAM File I/O
    /// \{

    /// \brief Helper method for access to underlying BGZF stream pointer.
    ///
    /// Useful for derived readers' contact points with htslib methods.
    ///
    /// \returns BGZF stream pointer
    ///
    BGZF* Bgzf() const;

    /// \brief Performs the actual raw read of the next record from the BAM
    ///        file.
    ///
    /// Default implementation will read records, sequentially, until EOF.
    /// Derived readers may use additional criteria to decide which record is
    ///  "next" and when reading is done.
    ///
    /// Return value should be equivalent to htslib's bam_read1():
    ///     >= 0 : normal
    ///       -1 : EOF (not an error)
    ///     < -1 : error
    ///
    /// \param[in]  bgzf BGZF stream pointer
    /// \param[out] b    %BAM record pointer
    /// \returns integer status code, see description
    ///
    virtual int ReadRawData(BGZF* bgzf, bam1_t* b);

    /// \}

private:
    std::unique_ptr<internal::BamReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BAMREADER_H
