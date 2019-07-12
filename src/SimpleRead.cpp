// File Description
/// \file SimpleRead.cpp
/// \brief Implements the SimpleRead class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SimpleRead.h"

#include <pbcopper/data/Clipping.h>

namespace PacBio {
namespace BAM {

void ClipToQuery(SimpleRead& read, Position start, Position end)
{
    PacBio::Data::ClipToQuery(read, start, end);
}

void ClipToQuery(MappedSimpleRead& read, Position start, Position end)
{
    PacBio::Data::ClipToQuery(read, start, end);
}

void ClipToReference(MappedSimpleRead& read, Position start, Position end,
                     bool exciseFlankingInserts)
{
    PacBio::Data::ClipToReference(read, start, end, exciseFlankingInserts);
}

}  // namespace BAM
}  // namespace PacBio
