// File Description
//
// Author: Derek Barnett, David Seifert

#ifndef PBBAM_MAKE_UNIQUE_H
#define PBBAM_MAKE_UNIQUE_H

// Only include if in C++11 mode or if using GCC 4.8
//
#if (__cplusplus <= 201103L) || ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8))

#include <cstddef>
#include <memory>

namespace std {

template <typename T, typename... Args>
inline std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}  // namespace std

#endif  // <= C++11

#endif  // PBBAM_MAKE_UNIQUE_H
