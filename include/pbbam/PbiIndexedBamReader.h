// File Description
/// \file PbiIndexedBamReader.h
/// \brief Defines the PbiIndexedBamReader class.
//
// Author: Derek Barnett

#ifndef PBIINDEXEDBAMREADER_H
#define PBIINDEXEDBAMREADER_H

#include <string>
#include "pbbam/BamFile.h"
#include "pbbam/BamReader.h"
#include "pbbam/PbiBasicTypes.h"
#include "pbbam/PbiFilter.h"

namespace PacBio {
namespace BAM {

namespace internal {
struct PbiIndexedBamReaderPrivate;
}

/// \brief The PbiIndexedBamReader class provides read-only iteration over %BAM
///        records, limited to some filtering criteria.
///
/// The PacBio BAM index (*.pbi) is used to allow random-access operations.
///
class PBBAM_EXPORT PbiIndexedBamReader : public BamReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Constructs %BAM reader, with an initial filter.
    ///
    /// All reads that satisfy the filter will be available.
    ///
    /// \param[in] filter       PbiFilter or compatible object
    /// \param[in] bamFilename  input %BAM filename
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(PbiFilter filter, const std::string& bamFilename);

    /// \brief Constructs %BAM reader, with an initial filter.
    ///
    /// All reads that satisfy the filter will be available.
    ///
    /// \param[in] filter       PbiFilter or compatible object
    /// \param[in] bamFile      input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(PbiFilter filter, BamFile bamFile);

    /// \brief Constructs %BAM reader, with no initial filter.
    ///
    /// Useful for delaying either specifying the filtering criteria or
    /// performing the PBI lookups.
    ///
    /// \param[in] bamFilename  input %BAM filename
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(const std::string& bamFilename);

    /// \brief Constructs %BAM reader, with no initial filter.
    ///
    /// Useful for delaying either specifying the filtering criteria or
    /// performing the PBI lookups.
    ///
    /// \param[in] bamFile      input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(BamFile bamFile);

    ~PbiIndexedBamReader() override;

    /// \}

public:
    /// \name Filtering & Index Data
    /// \{

    /// \returns the current filter active on this reader
    const PbiFilter& Filter() const;

    uint32_t NumReads() const;

    //    /// \returns the reader's underlying index data
    //    const PbiIndex& Index() const;

public:
    /// \brief Sets a new filter on the reader.
    ///
    /// \param[in] filter
    /// \returns reference to this reader
    ///
    PbiIndexedBamReader& Filter(PbiFilter filter);

    /// \}

protected:
    int ReadRawData(BGZF* bgzf, bam1_t* b) override;

private:
    std::unique_ptr<internal::PbiIndexedBamReaderPrivate> d_;
};

}  // namespace internal
}  // namespace BAM

#endif  // PBIINDEXEDBAMREADER_H
