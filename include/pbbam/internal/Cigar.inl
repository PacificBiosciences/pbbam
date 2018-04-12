// File Description
/// \file Cigar.inl
/// \brief Inline implemenations for the Cigar class.
//
// Author: Derek Barnett

#include "pbbam/Cigar.h"

namespace PacBio {
namespace BAM {

inline Cigar Cigar::FromStdString(const std::string& stdString)
{ return Cigar(stdString); }

} // namespace BAM
} // namespace PacBio
