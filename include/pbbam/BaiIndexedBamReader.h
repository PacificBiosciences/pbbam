// File Description
/// \file BaiIndexedBamReader.h
/// \brief Defines the BaiIndexedBamReader class.
//
// Author: Derek Barnett

#ifndef BAIINDEXEDBAMREADER_H
#define BAIINDEXEDBAMREADER_H

#include "pbbam/BamFile.h"
#include "pbbam/BamReader.h"
#include "pbbam/GenomicInterval.h"

namespace PacBio {
namespace BAM {

/// \brief The BaiIndexedBamReader class provides read-only iteration over %BAM
///        records, bounded by a particular genomic interval.
///
/// The SAM/BAM standard index (*.bai) is used to allow random-access operations.
///
class PBBAM_EXPORT BaiIndexedBamReader : public BamReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Constructs %BAM reader, that can be queried on genomic interval.
    ///
    /// \param filename input %BAM filename
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    explicit BaiIndexedBamReader(std::string filename);

    /// \brief Constructs %BAM reader, that can be queried on genomic interval.
    ///
    /// \param[in] bamFile   input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    explicit BaiIndexedBamReader(BamFile bamFile);

    /// \brief Constructs %BAM reader, bounded by a genomic interval.
    ///
    /// All reads that overlap the interval will be available.
    ///
    /// \param[in] interval iteration will be bounded by this GenomicInterval.
    /// \param[in] filename input %BAM filename
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    BaiIndexedBamReader(const GenomicInterval& interval, std::string filename);

    /// \brief Constructs %BAM reader, bounded by a genomic interval.
    ///
    /// All reads that overlap the interval will be available.
    ///
    /// \param[in] interval     iteration will be bounded by this GenomicInterval.
    /// \param[in] bamFile      input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    BaiIndexedBamReader(const GenomicInterval& interval, BamFile bamFile);

    /// \}

public:
    /// \name Random-Access
    /// \{

    /// \returns the underlying BamFile
    const BamFile& File() const;

    /// \returns the current GenomicInterval in use by this reader
    const GenomicInterval& Interval() const;

    /// \brief Sets a new genomic interval on the reader.
    ///
    /// \param[in] interval
    /// \returns reference to this reader
    ///
    BaiIndexedBamReader& Interval(const GenomicInterval& interval);

    /// \}

protected:
    int ReadRawData(BGZF* bgzf, bam1_t* b) override;

private:
    class BaiIndexedBamReaderPrivate;
    std::unique_ptr<BaiIndexedBamReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BAIINDEXEDBAMREADER_H
