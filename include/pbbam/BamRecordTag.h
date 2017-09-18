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
/// \file BamRecordTag.h
/// \brief Defines the BamRecordTag enum.
//
// Author: Derek Barnett

#ifndef BAMRECORDTAG_H
#define BAMRECORDTAG_H

namespace PacBio {
namespace BAM {

enum class BamRecordTag
{
    ALT_LABEL_QV
  , ALT_LABEL_TAG
  , BARCODE_QUALITY
  , BARCODES
  , CONTEXT_FLAGS
  , DELETION_QV
  , DELETION_TAG
  , HOLE_NUMBER
  , INSERTION_QV
  , IPD
  , LABEL_QV
  , MERGE_QV
  , NUM_PASSES
  , PKMEAN
  , PKMEAN_2
  , PKMID
  , PKMID_2
  , PRE_PULSE_FRAMES
  , PULSE_CALL
  , PULSE_CALL_WIDTH
  , PULSE_EXCLUSION
  , PULSE_MERGE_QV
  , PULSE_WIDTH
  , QUERY_END
  , QUERY_START
  , READ_ACCURACY
  , READ_GROUP
  , SCRAP_REGION_TYPE
  , SCRAP_ZMW_TYPE
  , SNR
  , START_FRAME
  , SUBSTITUTION_QV
  , SUBSTITUTION_TAG

  //
  // not tags per se, but faking these here to simplify data fetching
  //
  , QUAL
  , SEQ
};

} // namespace BAM
} // namespace PacBio

#endif // BAMRECORDTAG_H
