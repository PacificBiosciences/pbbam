#ifndef PBBAM_FASTQSEQUENCE_H
#define PBBAM_FASTQSEQUENCE_H

#include <pbbam/Config.h>

#include <pbbam/FastaSequence.h>

#include <pbcopper/data/QualityValues.h>

#include <string>

namespace PacBio {
namespace BAM {

///
/// \brief The FastqSequence class represents a FASTQ record (name, bases, and
///        qualities)
///
class FastqSequence : public FastaSequence
{
public:
    /// \name Constructors & Related Methods
    /// \{

    ///
    /// \brief FastaSequence
    /// \param name
    /// \param bases
    /// \param qualities
    ///
    explicit FastqSequence(std::string name, std::string bases, Data::QualityValues qualities);

    ///
    /// \brief FastaSequence
    /// \param name
    /// \param bases
    /// \param qualities
    ///
    explicit FastqSequence(std::string name, std::string bases, std::string qualities);

    FastqSequence() = default;

    /// \}

public:
    /// \name Attributes
    /// \{

    ///
    /// \brief Qualities
    /// \return
    ///
    const Data::QualityValues& Qualities() const;

    ///
    /// \brief
    ///
    /// \param quals
    /// \return FastqSequence
    ///
    FastqSequence& Qualities(Data::QualityValues quals);

    /// \}

    ///
    /// Return average base quality.
    ///
    /// \throw std::runtime_error on empty sequence
    ///
    float AverageBaseQuality() const;

    bool operator==(const FastqSequence& other) const noexcept;
    bool operator!=(const FastqSequence& other) const noexcept;

private:
    Data::QualityValues qualities_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FASTQSEQUENCE_H
