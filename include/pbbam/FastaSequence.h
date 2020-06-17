// File Description
/// \file FastaSequence.h
/// \brief Defines the FastaSequence class.
//
// Author: Derek Barnett

#ifndef FASTASEQUENCE_H
#define FASTASEQUENCE_H

#include "pbbam/Config.h"

#include <string>

namespace PacBio {
namespace BAM {

///
/// \brief The FastaSequence class represents a FASTA record (name & bases)
///
class FastaSequence
{
public:
    /// \name Constructors & Related Methods
    /// \{

    ///
    /// \brief FastaSequence
    /// \param name
    /// \param bases
    ///
    explicit FastaSequence(std::string name, std::string bases);

    FastaSequence() = default;

    /// \}

public:
    /// \name Attributes
    /// \{

    ///
    /// \brief Name
    /// \return
    ///
    const std::string& Name() const;

    ///
    /// \brief
    ///
    /// \param name
    /// \return FastaSequence&
    ///
    FastaSequence& Name(std::string name);

    ///
    /// \brief Bases
    /// \return
    ///
    const std::string& Bases() const;

    ///
    /// \brief
    ///
    /// \param bases
    /// \return FastaSequence&
    ///
    FastaSequence& Bases(std::string bases);

    /// \}

    bool operator==(const FastaSequence& other) const noexcept;
    bool operator!=(const FastaSequence& other) const noexcept;

private:
    std::string name_;
    std::string bases_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTASEQUENCE_H
