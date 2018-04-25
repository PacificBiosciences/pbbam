// File Description
/// \file Orientation.h
/// \brief Defines the Orientation enum.
//
// Author: Derek Barnett

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief This enum defines the orientations recognized by BamRecord, for
///        presenting "per-base" data.
///
/// Orientation::NATIVE indicates that data should be presented in the subread's
/// original form.
///
/// Orientation::GENOMIC indicates that data should be presented relative to
/// genomic forward strand. This means that data will be reversed (or
/// reverse-complemented) if the subread was aligned to the reverse strand.
///
enum class Orientation
{
    NATIVE,  ///< Present data in 'raw' original orientation, regardless of aligned Strand
    GENOMIC  ///< Present data in aligned orientation, always relative to Strand::FORWARD.
};

}  // namespace BAM
}  // namespace PacBio

#endif  // ORIENTATION_H
