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

using PacBio::Data::ADAPTER_AFTER;
using PacBio::Data::ADAPTER_AFTER_BAD;
using PacBio::Data::ADAPTER_BEFORE;
using PacBio::Data::ADAPTER_BEFORE_BAD;
using PacBio::Data::BARCODE_AFTER;
using PacBio::Data::BARCODE_BEFORE;
using PacBio::Data::FORWARD_PASS;
using PacBio::Data::NO_LOCAL_CONTEXT;
using PacBio::Data::REVERSE_PASS;

}  // namespace BAM
}  // namespace PacBio

#endif  // LOCALCONTEXTFLAGS_H
