// Author: Derek Barnett

#ifndef MOVEAPPEND_H
#define MOVEAPPEND_H

#include <iterator>
#include <utility>
#include <vector>

namespace PacBio {
namespace BAM {

// \brief Appends content of src vector to dst vector using move semantics.
///
/// \param[in]     src  Input vector that will be empty after execution
/// \param[in,out] dst  Output vector that will be appended to
///
template <typename T>
inline void MoveAppend(std::vector<T>& src, std::vector<T>& dst) noexcept
{
    if (dst.empty()) {
        dst = std::move(src);
    } else {
        dst.reserve(dst.size() + src.size());
        std::move(src.begin(), src.end(), std::back_inserter(dst));
        src.clear();
    }
}

/// \brief Appends content of src vector to dst vector using move semantics.
///
/// \param[in]     src  Input vector via perfect forwarding
/// \param[in,out] dst  Output vector that will be appended to
///
template <typename T>
inline void MoveAppend(std::vector<T>&& src, std::vector<T>& dst) noexcept
{
    if (dst.empty()) {
        dst = std::move(src);
    } else {
        dst.reserve(dst.size() + src.size());
        std::move(src.begin(), src.end(), std::back_inserter(dst));
        src.clear();
    }
}

}  // namespace BAM
}  // namespace PacBio

#endif  // MOVEAPPEND_H
