// File Description
/// \file Strand.h
/// \brief Defines the Strand enum.
//
// Author: Derek Barnett

#ifndef STRAND_H
#define STRAND_H

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief This enum defines the strand orientations used for reporting
///        alignment-related information.
///
enum class Strand
{
    FORWARD,  ///< Forward strand
    REVERSE   ///< Reverse strand
};

}  // namespace BAM
}  // namespace PacBio

#endif  // STRAND_H
