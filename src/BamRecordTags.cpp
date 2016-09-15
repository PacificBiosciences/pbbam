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
/// \file BamRecordTags.h
/// \brief Implements the BamRecordTags utility class.
//
// Author: Derek Barnett

#include "BamRecordTags.h"
#include "EnumClassHash.h"
#include <unordered_map>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

const BamRecordTags::TagLookupType BamRecordTags::tagLookup =
{
    //     enum name                   label  isPulse?
    //     ---------                   -----  --------
    { BamRecordTag::ALT_LABEL_QV,      {"pv", true}  },
    { BamRecordTag::ALT_LABEL_TAG,     {"pt", true}  },
    { BamRecordTag::BARCODE_QUALITY,   {"bq", false} },
    { BamRecordTag::BARCODES,          {"bc", false} },
    { BamRecordTag::CONTEXT_FLAGS,     {"cx", false} },
    { BamRecordTag::DELETION_QV,       {"dq", false} },
    { BamRecordTag::DELETION_TAG,      {"dt", false} },
    { BamRecordTag::HOLE_NUMBER,       {"zm", false} },
    { BamRecordTag::INSERTION_QV,      {"iq", false} },
    { BamRecordTag::IPD,               {"ip", false} },
    { BamRecordTag::LABEL_QV,          {"pq", true}  },
    { BamRecordTag::MERGE_QV,          {"mq", false} },
    { BamRecordTag::NUM_PASSES,        {"np", false} },
    { BamRecordTag::PKMEAN,            {"pa", true}  },
    { BamRecordTag::PKMEAN_2,          {"ps", true}  },
    { BamRecordTag::PKMID,             {"pm", true}  },
    { BamRecordTag::PKMID_2,           {"pi", true}  },
    { BamRecordTag::PRE_PULSE_FRAMES,  {"pd", true}  },
    { BamRecordTag::PULSE_CALL,        {"pc", true}  },
    { BamRecordTag::PULSE_CALL_WIDTH,  {"px", true}  },
    { BamRecordTag::PULSE_MERGE_QV,    {"pg", true}  },
    { BamRecordTag::PULSE_WIDTH,       {"pw", false} }, // 'pulse' in the name; but stored per-base, not per-pulse
    { BamRecordTag::QUERY_END,         {"qe", false} },
    { BamRecordTag::QUERY_START,       {"qs", false} },
    { BamRecordTag::READ_ACCURACY,     {"rq", false} },
    { BamRecordTag::READ_GROUP,        {"RG", false} },
    { BamRecordTag::SCRAP_REGION_TYPE, {"sc", false} },
    { BamRecordTag::SCRAP_ZMW_TYPE,    {"sz", false} },
    { BamRecordTag::SNR,               {"sn", false} },
    { BamRecordTag::START_FRAME,       {"sf", true}  },
    { BamRecordTag::SUBSTITUTION_QV,   {"sq", false} },
    { BamRecordTag::SUBSTITUTION_TAG,  {"st", false} },

    // faux tags
    { BamRecordTag::SEQ,  {"  ",  false} },
    { BamRecordTag::QUAL, {"  ", false} }
};

} // namespace internal
} // namespace BAM
} // namespace PacBio
