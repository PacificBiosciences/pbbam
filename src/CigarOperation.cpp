// File Description
/// \file CigarOperation.cpp
/// \brief Implements the CigarOperation class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/CigarOperation.h"

#include <htslib/sam.h>

namespace PacBio {
namespace BAM {

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

char CigarOperation::TypeToChar(const CigarOperationType type)
{
    return bam_cigar_opchr(static_cast<int>(type));
}

bool CigarOperation::validate_ = true;

}  // namespace BAM
}  // namespace PacBio
