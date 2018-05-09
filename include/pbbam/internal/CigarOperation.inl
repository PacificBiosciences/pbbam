// File Description
/// \file CigarOperation.inl
/// \brief Inline implemenations for the CigarOperation class.
//
// Author: Derek Barnett

#include "pbbam/CigarOperation.h"

namespace PacBio {
namespace BAM {

inline CigarOperation::CigarOperation(char c, uint32_t length)
    : type_{CigarOperation::CharToType(c)}
    , length_{length}
{
    #ifndef PBBAM_PERMISSIVE_CIGAR
        if (validate_ && (type_ == CigarOperationType::ALIGNMENT_MATCH))
            throw std::runtime_error{"CIGAR operation 'M' is not allowed in PacBio BAM files. Use 'X/=' instead."};
    #endif
}

inline CigarOperation::CigarOperation(CigarOperationType op, uint32_t length)
    : type_{op}
    , length_{length}
{
    #ifndef PBBAM_PERMISSIVE_CIGAR
        if (validate_ && (type_ == CigarOperationType::ALIGNMENT_MATCH))
            throw std::runtime_error{"CIGAR operation 'M' is not allowed in PacBio BAM files. Use 'X/=' instead."};
    #endif
}

inline uint32_t CigarOperation::Length() const
{ return length_; }

inline CigarOperation& CigarOperation::Length(const uint32_t length)
{ length_ = length; return *this; }

inline CigarOperationType CigarOperation::Type() const
{ return type_; }

inline CigarOperation &CigarOperation::Type(const CigarOperationType opType)
{ type_ = opType; return *this; }

inline char CigarOperation::Char() const
{ return CigarOperation::TypeToChar(type_); }

inline CigarOperation &CigarOperation::Char(const char opChar)
{ type_ = CigarOperation::CharToType(opChar); return *this; }

inline bool CigarOperation::operator==(const CigarOperation& other) const
{ return type_ == other.type_ && length_ == other.length_; }

inline bool CigarOperation::operator!=(const CigarOperation& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio
