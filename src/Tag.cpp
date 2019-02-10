// File Description
/// \file Tag.cpp
/// \brief Implements the Tag class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Tag.h"

namespace PacBio {
namespace BAM {

Tag::Tag() = default;

Tag::Tag(const Tag&) = default;

Tag::Tag(Tag&&) = default;

Tag& Tag::operator=(const Tag&) = default;

Tag& Tag::operator=(Tag&&) = default;

Tag::~Tag() = default;

}  // namespace BAM
}  // namespace PacBio