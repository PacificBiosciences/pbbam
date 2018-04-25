// File Description
/// \file LocalContextFlags.h
/// \brief Defines the LocalContextFlags enum & helper method(s).
//
// Author: Lance Hepler

#ifndef LOCALCONTEXTFLAGS_H
#define LOCALCONTEXTFLAGS_H

#include <cstdint>

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief The LocalContextFlags enum defines the flags that can be used
///        to describe a subread's "local context", i.e. whether it is
///        flanked by barcodes/adapters or its pass orientation.
///
enum LocalContextFlags : uint8_t
{
    NO_LOCAL_CONTEXT = 0,  ///< No context information available
    ADAPTER_BEFORE = 1,    ///< Adapter precedes subread
    ADAPTER_AFTER = 2,     ///< Adapter follows subread
    BARCODE_BEFORE = 4,    ///< Barcode precedes subread
    BARCODE_AFTER = 8,     ///< Barcode follows subread
    FORWARD_PASS = 16,     ///< Subread's orientation is 'forward pass'
    REVERSE_PASS = 32      ///< Subread's orientation is 'reverse pass'
};

/// \returns a LocalContextFlags value containing the result of the bitwise-OR
///          operation of \p lhs and \p rhs.
// constexpr is implicitly inline
constexpr LocalContextFlags operator|(const LocalContextFlags lhs, const LocalContextFlags rhs)
{
    return static_cast<LocalContextFlags>(static_cast<int>(lhs) | static_cast<int>(rhs));
}

}  // namespace BAM
}  // namespace PacBio

#endif  // LOCALCONTEXTFLAGS_H
