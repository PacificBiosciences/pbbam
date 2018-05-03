// File Description
/// \file QualityValue.h
/// \brief Defines the QualityValue class.
//
// Author: Derek Barnett

#ifndef QUALITYVALUE_H
#define QUALITYVALUE_H

#include <cstdint>
#include <string>
#include <vector>
#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief The QualityValue class represents a FASTQ-compatible quality value.
///
/// Integers are clamped to [0, 93] (corresponding to ASCII printable chars
/// [!-~]).
///
/// Use QualityValue::FromFastq for constructing entries from FASTQ encoding
/// characters. Otherwise, the resulting QualityValue will be interpreted using
/// the character's numeric value (ignoring the FASTQ offset of 33).
///
class PBBAM_EXPORT QualityValue
{
public:
    static const uint8_t MAX;

public:
    /// \name Conversion Methods
    /// \{

    /// \brief Creates a QualityValue from a FASTQ-encoding character.
    ///
    /// \param[in] c    FASTQ character
    /// \returns quality value representing (c - 33)
    ///
    static QualityValue FromFastq(const char c);

    /// \}

public:
    /// \name Constructors & Related Methods
    ///  \{

    /// \brief Creates a QualityValue with specified value.
    ///
    /// \param[in] value    quality value
    ///
    QualityValue(const uint8_t value = 0);

    QualityValue(const QualityValue&) = default;
    QualityValue(QualityValue&&) = default;
    QualityValue& operator=(const QualityValue&) = default;
    QualityValue& operator=(QualityValue&&) = default;
    ~QualityValue() = default;

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// \returns the FASTQ-encoding char for this QualityValue
    char Fastq() const;

    /// \returns the integer value of this QualityValue
    operator uint8_t() const;

    /// \}

private:
    uint8_t value_;
};

}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/QualityValue.inl"

#endif  // QUALITYVALUE_H
