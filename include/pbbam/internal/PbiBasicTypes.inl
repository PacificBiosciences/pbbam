#ifndef PBBAM_PBIBASICTYPES_INL
#define PBBAM_PBIBASICTYPES_INL

#include <pbbam/Config.h>

#include <pbbam/PbiBasicTypes.h>

#include <tuple>

namespace PacBio {
namespace BAM {

inline IndexResultBlock::IndexResultBlock(size_t idx, size_t numReads)
    : firstIndex_{idx}, numReads_{numReads}
{
}

inline bool IndexResultBlock::operator==(const IndexResultBlock& other) const noexcept
{
    return std::tie(firstIndex_, numReads_, virtualOffset_) ==
           std::tie(other.firstIndex_, other.numReads_, other.virtualOffset_);
}

inline bool IndexResultBlock::operator!=(const IndexResultBlock& other) const noexcept
{
    return !(*this == other);
}

}  // namespace BAM
}  // namespace PacBio

#endif // PBBAM_PBIBASICTYPES_INL
