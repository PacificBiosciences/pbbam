// File Description
/// \file SImpleRead.h
/// \brief Defines the SimpleRead class.
//
// Author: Derek Barnett

#ifndef SIMPLEREAD_H
#define SIMPLEREAD_H

#include "pbbam/Config.h"

// #include <string>

// #include <boost/optional.hpp>

// #include "pbbam/BamRecord.h"
// #include "pbbam/Cigar.h"
// #include "pbbam/Frames.h"
// #include "pbbam/Position.h"
// #include "pbbam/QualityValues.h"
// #include "pbbam/SNR.h"
// #include "pbbam/Strand.h"

#include <pbcopper/data/MappedSimpleRead.h>
#include <pbcopper/data/SimpleRead.h>

#include "pbbam/Position.h"

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
void ClipToQuery(SimpleRead& read, Position start, Position end);

///
/// \brief
///
/// \param read
/// \param start
/// \param end
///
void ClipToQuery(MappedSimpleRead& read, Position start, Position end);

///
/// \brief
///
/// \param read
/// \param start
/// \param end
/// \param exciseFlankingInserts
///
void ClipToReference(MappedSimpleRead& read, Position start, Position end,
                     bool exciseFlankingInserts);

}  // namespace BAM
}  // namespace PacBio

#endif  // SIMPLEREAD_H
