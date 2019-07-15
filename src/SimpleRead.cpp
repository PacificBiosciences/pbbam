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

void ClipToQuery(Data::SimpleRead& read, Data::Position start, Data::Position end)
{
    Data::ClipToQuery(read, start, end);
}

void ClipToQuery(Data::MappedSimpleRead& read, Data::Position start, Data::Position end)
{
    Data::ClipToQuery(read, start, end);
}

void ClipToReference(Data::MappedSimpleRead& read, Data::Position start, Data::Position end,
                     bool exciseFlankingInserts)
{
    Data::ClipToReference(read, start, end, exciseFlankingInserts);
}

}  // namespace BAM
}  // namespace PacBio
