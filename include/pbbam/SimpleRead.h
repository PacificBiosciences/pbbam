// File Description
/// \file SImpleRead.h
/// \brief Defines the SimpleRead class.
//
// Author: Derek Barnett

#ifndef SIMPLEREAD_H
#define SIMPLEREAD_H

#include "pbbam/Config.h"

#include <pbcopper/data/MappedSimpleRead.h>
#include <pbcopper/data/Position.h>
#include <pbcopper/data/SimpleRead.h>

#ifndef PBBAM_NODEPRECATED_API

namespace PacBio {
namespace BAM {

using SimpleRead = PacBio::Data::SimpleRead;
using MappedSimpleRead = PacBio::Data::MappedSimpleRead;

///
/// \brief
///
/// \param read
/// \param start
/// \param end
///
void ClipToQuery(Data::SimpleRead& read, Data::Position start, Data::Position end);

///
/// \brief
///
/// \param read
/// \param start
/// \param end
///
void ClipToQuery(Data::MappedSimpleRead& read, Data::Position start, Data::Position end);

///
/// \brief
///
/// \param read
/// \param start
/// \param end
/// \param exciseFlankingInserts
///
void ClipToReference(Data::MappedSimpleRead& read, Data::Position start, Data::Position end,
                     bool exciseFlankingInserts);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_NODEPRECATED_API

#endif  // SIMPLEREAD_H
