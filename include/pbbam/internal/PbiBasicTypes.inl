// File Description
/// \file PbiBasicTypes.inl
/// \brief Inline implementations for the basic data structures used in PBI lookups.
//
// Author: Derek Barnett

#include "pbbam/PbiBasicTypes.h"

namespace PacBio {
namespace BAM {

inline IndexResultBlock::IndexResultBlock(size_t idx, size_t numReads)
    : firstIndex_{idx}
    , numReads_{numReads}
{ }

inline bool IndexResultBlock::operator==(const IndexResultBlock& other) const
{
    return firstIndex_ == other.firstIndex_ &&
           numReads_ == other.numReads_ &&
           virtualOffset_ == other.virtualOffset_;
}

inline bool IndexResultBlock::operator!=(const IndexResultBlock& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio
