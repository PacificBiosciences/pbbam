// Copyright (c) 2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// File Description
/// \file RecordType.h
/// \brief Defines the RecordType enum.
//
// Author: Derek Barnett

#ifndef RECORDTYPE_H
#define RECORDTYPE_H

namespace PacBio {
namespace BAM {

/// \brief This enum defines the possible PacBio BAM record types.
///
/// \sa ReadGroupInfo::ReadType
///
enum class RecordType
{
    ZMW         ///< Polymerase read
  , HQREGION    ///< High-quality region
  , SUBREAD     ///< Subread (
  , CCS         ///< Circular consensus sequence
  , SCRAP       ///< Additional sequence (barcodes, adapters, etc.)
  , UNKNOWN     ///< Unknown read type

  , POLYMERASE = ZMW ///< \deprecated as of PacBio BAM spec v 3.0.4 (use RecordType::ZMW instead)
};

} // namespace BAM
} // namespace PacBio

#endif // RECORDTYPE_H
