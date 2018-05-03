// File Description
/// \file AlignmentPrinter.h
/// \brief Defines the AlignmentPrinter class.
//
// Author: Armin Töpfer

#ifndef ALIGNMENTPRINTER_H
#define ALIGNMENTPRINTER_H

#include <memory>
#include <string>

#include "pbbam/BamRecord.h"
#include "pbbam/IndexedFastaReader.h"
#include "pbbam/Orientation.h"

namespace PacBio {
namespace BAM {

class BamRecord;

/// \brief The AlignmentPrinter class "pretty-prints" an alignment with respect
///        to its associated reference sequence.
///
/// Example output:
/// \verbinclude plaintext/AlignmentPrinterOutput.txt
///
class AlignmentPrinter
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// Constructs the alignment printer with an associated FASTA file reader.
    ///
    /// \param[in] ifr FASTA reader
    ///
    /// \throws std::runtime_error if FASTA file cannot be opened for reading.
    ///
    AlignmentPrinter(const IndexedFastaReader& ifr);

    AlignmentPrinter() = delete;
    AlignmentPrinter(const AlignmentPrinter&) = delete;
    AlignmentPrinter(AlignmentPrinter&&) = default;
    AlignmentPrinter& operator=(const AlignmentPrinter&) = delete;
    AlignmentPrinter& operator=(AlignmentPrinter&&) = default;
    ~AlignmentPrinter() = default;

    /// \}

public:
    /// \name Printing
    /// \{

    /// Pretty-prints an aligned BamRecord to std::string.
    ///
    /// \note The current implementation includes ANSI escape sequences for
    ///       coloring terminal output. Future versions of this method will
    ///       likely make this optional.
    ///
    /// \returns formatted string containing the alignment and summary
    ///          information
    ///
    std::string Print(const BamRecord& record,
                      const Orientation orientation = Orientation::GENOMIC);

    /// \}

private:
    const std::unique_ptr<IndexedFastaReader> ifr_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // ALIGNMENTPRINTER_H
