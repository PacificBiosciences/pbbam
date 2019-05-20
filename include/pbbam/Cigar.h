// File Description
/// \file Cigar.h
/// \brief Defines the Cigar class.
//
// Author: Derek Barnett

#ifndef CIGAR_H
#define CIGAR_H

#include "pbbam/Config.h"

#include <string>
#include <vector>

#include "pbbam/CigarOperation.h"

namespace PacBio {
namespace BAM {

/// \brief The Cigar class represents the CIGAR string used to report alignment
///        charateristics in SAM/BAM.
///
/// \note Use of the 'M' operator is forbidden in PacBio BAMs. See
///       CigarOperationType description for more information.
///
/// \sa https://samtools.github.io/hts-specs/SAMv1.pdf for more information on CIGAR in general.
///
class PBBAM_EXPORT Cigar : public std::vector<CigarOperation>
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a Cigar object from SAM/BAM string input
    ///
    /// \param [in] stdString   SAM/BAM formatted CIGAR data
    /// \returns a Cigar object representing the input data
    ///
    /// \note This class may be removed from the public API in the future,
    ///       as the constructor taking a std::string accomplishes the same end.
    ///
    static Cigar FromStdString(const std::string& stdString);

    /// \brief Creates an empty Cigar.
    Cigar();

    /// \brief Creates a Cigar object from SAM/BAM string input
    ///
    /// \param [in] cigarString   SAM/BAM formatted CIGAR data
    ///
    Cigar(const std::string& cigarString);

    Cigar(const Cigar&);
    Cigar(Cigar&&) noexcept;
    Cigar& operator=(const Cigar&);
    Cigar& operator=(Cigar&&) noexcept;
    ~Cigar();

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// Converts Cigar object data to SAM/BAM formatted string
    ///
    /// \returns SAM/BAM formatted std::string
    ///
    std::string ToStdString() const;

    /// \}
};

///
/// \brief
///
/// \param cigar
/// \return size_t
///
size_t ReferenceLength(const Cigar& cigar);

}  // namespace BAM
}  // namespace PacBio

#endif  // CIGAR_H
