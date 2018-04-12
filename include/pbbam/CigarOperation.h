// File Description
/// \file CigarOperation.h
/// \brief Defines the CigarOperationType enum & CigarOperation class.
//
// Author: Derek Barnett

#ifndef CIGAROPERATION_H
#define CIGAROPERATION_H

#include <cstdint>
#include <stdexcept>
#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

namespace pbbamify {
class Settings;
}

/// \brief Describes a CIGAR operation.
///
/// Bracketed character is the corresponding SAM/BAM character code.
///
/// \warning ALIGNMENT_MATCH ('M') is included in this enum to maintain
///          consistency with htslib. However, as of PacBio BAM spec version
///          3.0b7, this CIGAR operation \b forbidden. Any attempt to read or
///          write a record containing this operation will trigger a
///          std::runtime_error. SEQUENCE_MATCH('=) or SEQUENCE_MISMATCH('X')
///          should be used instead.
///
enum class CigarOperationType
{
    UNKNOWN_OP = -1,      ///< unknown/invalid CIGAR operator
    ALIGNMENT_MATCH = 0,  ///< alignment match (can be a sequence match or mismatch) [M]
    INSERTION,            ///< insertion to the reference [I]
    DELETION,             ///< deletion from the reference [D]
    REFERENCE_SKIP,       ///< skipped region from the reference [N]
    SOFT_CLIP,            ///< soft clipping (clipped sequences present in SEQ) [S]
    HARD_CLIP = 5,        ///< hard clipping (clipped sequences NOT present in SEQ) [H]
    PADDING,              ///< padding (silent deletion from padded reference) [P]
    SEQUENCE_MATCH,       ///< sequence match [=]
    SEQUENCE_MISMATCH     ///< sequence mismatch [X]
};

/// \brief The CigarOperation class represents a single CIGAR operation
///        (consisting of a type & length).
///
class PBBAM_EXPORT CigarOperation
{
public:
    /// \name Operation Type Conversion Methods
    /// \{

    /// Convert between CigarOperationType enum & SAM/BAM character code.
    ///
    /// \param[in] type CigarOperationType value
    /// \returns SAM/BAM character code
    static char TypeToChar(const CigarOperationType type);

    /// Convert between CigarOperationType enum & SAM/BAM character code.
    ///
    /// \param[in] c SAM/BAM character code
    /// \returns CigarOperationType value
    static CigarOperationType CharToType(const char c);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    CigarOperation() = default;
    CigarOperation(char c, uint32_t length);
    CigarOperation(CigarOperationType op, uint32_t length);

    CigarOperation(const CigarOperation&) = default;
    CigarOperation(CigarOperation&&) = default;
    CigarOperation& operator=(const CigarOperation&) = default;
    CigarOperation& operator=(CigarOperation&&) = default;
    ~CigarOperation() = default;

    /// \}

public:
    /// \returns operation type as SAM/BAM char code
    inline char Char() const;

    /// \returns operation length
    inline uint32_t Length() const;

    /// \returns operation type as CigarOperationType enum value
    inline CigarOperationType Type() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// Sets this operation type.
    ///
    /// \param[in] opChar SAM/BAM character code
    /// \returns reference to this operation
    inline CigarOperation& Char(const char opChar);

    /// Sets this operation length.
    ///
    /// \param[in] length
    /// \returns reference to this operation
    inline CigarOperation& Length(const uint32_t length);

    /// Sets this operation type.
    ///
    /// \param[in] opType CigarOperationType value
    /// \returns reference to this operation
    inline CigarOperation& Type(const CigarOperationType opType);

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    /// \returns true if both CIGAR operation type & length match
    inline bool operator==(const CigarOperation& other) const;

    /// \returns true if either CIGAR operation type or length differ
    inline bool operator!=(const CigarOperation& other) const;

    /// \}

private:
    CigarOperationType type_ = CigarOperationType::UNKNOWN_OP;
    uint32_t length_ = 0;

    // runtime disabling of Cigar validation
    static bool validate_;
    friend class pbbamify::Settings;
};

}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/CigarOperation.inl"

#endif  // CIGAROPERATION_H
