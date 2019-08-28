// File Description
/// \file LocalContextFlags.h
/// \brief Defines the LocalContextFlags enum & helper method(s).
//
// Author: Lance Hepler

#ifndef LOCALCONTEXTFLAGS_H
#define LOCALCONTEXTFLAGS_H

#include "pbbam/Config.h"

#include <pbcopper/data/LocalContextFlags.h>

namespace PacBio {
namespace BAM {

using LocalContextFlags PBBAM_DEPRECATED = PacBio::Data::LocalContextFlags;

// because LocalContextFlags was a C enum and not a
// C++11 enum class, we need to import all enumerations
// into the containing scope
using PacBio::Data::LocalContextFlags::NO_LOCAL_CONTEXT;
using PacBio::Data::LocalContextFlags::ADAPTER_BEFORE;
using PacBio::Data::LocalContextFlags::ADAPTER_AFTER;
using PacBio::Data::LocalContextFlags::BARCODE_BEFORE;
using PacBio::Data::LocalContextFlags::BARCODE_AFTER;
using PacBio::Data::LocalContextFlags::FORWARD_PASS;
using PacBio::Data::LocalContextFlags::REVERSE_PASS;
using PacBio::Data::LocalContextFlags::ADAPTER_BEFORE_BAD;
using PacBio::Data::LocalContextFlags::ADAPTER_AFTER_BAD;

}  // namespace BAM
}  // namespace PacBio

#endif  // LOCALCONTEXTFLAGS_H
