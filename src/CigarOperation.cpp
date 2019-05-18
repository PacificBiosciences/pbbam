// File Description
/// \file CigarOperation.cpp
/// \brief Implements the CigarOperation class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/CigarOperation.h"

#include <tuple>

#include <htslib/sam.h>

namespace PacBio {
namespace BAM {

bool CigarOperation::validate_ = true;

CigarOperation::CigarOperation() = default;

CigarOperation::CigarOperation(char c, uint32_t length)
    : type_{CigarOperation::CharToType(c)}, length_{length}
{
#ifndef PBBAM_PERMISSIVE_CIGAR
    if (validate_ && (type_ == CigarOperationType::ALIGNMENT_MATCH))
        throw std::runtime_error{
            "CIGAR operation 'M' is not allowed in PacBio BAM files. Use 'X/=' instead."};
#endif
}

CigarOperation::CigarOperation(CigarOperationType op, uint32_t length) : type_{op}, length_{length}
{
#ifndef PBBAM_PERMISSIVE_CIGAR
    if (validate_ && (type_ == CigarOperationType::ALIGNMENT_MATCH))
        throw std::runtime_error{
            "CIGAR operation 'M' is not allowed in PacBio BAM files. Use 'X/=' instead."};
#endif
}

CigarOperation::CigarOperation(const CigarOperation&) = default;

CigarOperation::CigarOperation(CigarOperation&&) noexcept = default;

CigarOperation& CigarOperation::operator=(const CigarOperation&) = default;

CigarOperation& CigarOperation::operator=(CigarOperation&&) noexcept = default;

CigarOperation::~CigarOperation() = default;

bool CigarOperation::operator==(const CigarOperation& other) const
{
    return std::tie(type_, length_) == std::tie(other.type_, other.length_);
}

bool CigarOperation::operator!=(const CigarOperation& other) const { return !(*this == other); }

char CigarOperation::Char() const { return CigarOperation::TypeToChar(type_); }

CigarOperation& CigarOperation::Char(const char opChar)
{
    type_ = CigarOperation::CharToType(opChar);
    return *this;
}

CigarOperationType CigarOperation::CharToType(const char c)
{
    switch (c) {
        case 'S':
            return CigarOperationType::SOFT_CLIP;
        case '=':
            return CigarOperationType::SEQUENCE_MATCH;
        case 'X':
            return CigarOperationType::SEQUENCE_MISMATCH;
        case 'I':
            return CigarOperationType::INSERTION;
        case 'D':
            return CigarOperationType::DELETION;
        case 'N':
            return CigarOperationType::REFERENCE_SKIP;
        case 'H':
            return CigarOperationType::HARD_CLIP;
        case 'P':
            return CigarOperationType::PADDING;
        case 'M':
            return CigarOperationType::ALIGNMENT_MATCH;
        default:
            return CigarOperationType::UNKNOWN_OP;
    }
}

uint32_t CigarOperation::Length() const { return length_; }

CigarOperation& CigarOperation::Length(const uint32_t length)
{
    length_ = length;
    return *this;
}

CigarOperationType CigarOperation::Type() const { return type_; }

CigarOperation& CigarOperation::Type(const CigarOperationType opType)
{
    type_ = opType;
    return *this;
}

char CigarOperation::TypeToChar(const CigarOperationType type)
{
    return bam_cigar_opchr(static_cast<int>(type));
}

}  // namespace BAM
}  // namespace PacBio
