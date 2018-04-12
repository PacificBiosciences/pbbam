// File Description
/// \file BamRecord.h
/// \brief Defines the BamRecord class.
//
// Author: Derek Barnett

#ifndef PBBAM_MAKE_UNIQUE_H
#define PBBAM_MAKE_UNIQUE_H

// Only include if in pre-C++14 mode
//
#if __cplusplus <= 201103L

#include <cstddef>
#include <memory>

namespace std {

template <typename T, typename... Args>
inline std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}  // namespace std

#endif  // < C++14

#endif  // PBBAM_MAKE_UNIQUE_H
