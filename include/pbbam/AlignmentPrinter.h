#ifndef PBBAM_ALIGNMENTPRINTER_H
#define PBBAM_ALIGNMENTPRINTER_H

#include <pbbam/Config.h>

#include <memory>
#include <string>

#include <pbcopper/data/Orientation.h>

#include <pbbam/BamRecord.h>
#include <pbbam/IndexedFastaReader.h>

namespace PacBio {
namespace BAM {

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
                      const Data::Orientation orientation = Data::Orientation::GENOMIC);

    /// \}

private:
    std::unique_ptr<IndexedFastaReader> ifr_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ALIGNMENTPRINTER_H
