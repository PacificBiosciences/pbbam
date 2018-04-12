// File Description
/// \file TagCollection.cpp
/// \brief Implements the TagCollection class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/TagCollection.h"

namespace PacBio {
namespace BAM {

bool TagCollection::Contains(const std::string& name) const { return count(name) != 0; }

}  // namespace BAM
}  // namespace PacBio
