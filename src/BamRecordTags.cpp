#include "PbbamInternalConfig.h"

#include "BamRecordTags.h"

#include <unordered_map>

namespace PacBio {
namespace BAM {

// clang-format off
const BamRecordTags::TagLookupType BamRecordTags::tagLookup =
{
    //     enum name                          label  isPulse?
    //     ---------                          -----  --------
    { BamRecordTag::ALT_LABEL_QV,             {"pv", true}  },
    { BamRecordTag::ALT_LABEL_TAG,            {"pt", true}  },
    { BamRecordTag::BARCODE_QUALITY,          {"bq", false} },
    { BamRecordTag::BARCODES,                 {"bc", false} },
    { BamRecordTag::CONTEXT_FLAGS,            {"cx", false} },
    { BamRecordTag::DELETION_QV,              {"dq", false} },
    { BamRecordTag::DELETION_TAG,             {"dt", false} },
    { BamRecordTag::FORWARD_IPD,              {"fi", false} },
    { BamRecordTag::FORWARD_PW,               {"fp", false} },
    { BamRecordTag::HOLE_NUMBER,              {"zm", false} },
    { BamRecordTag::INSERTION_QV,             {"iq", false} },
    { BamRecordTag::IPD,                      {"ip", false} },
    { BamRecordTag::LABEL_QV,                 {"pq", true}  },
    { BamRecordTag::LONG_CIGAR,               {"CG", false} },
    { BamRecordTag::MERGE_QV,                 {"mq", false} },
    { BamRecordTag::NUM_PASSES,               {"np", false} },
    { BamRecordTag::PKMEAN,                   {"pa", true}  },
    { BamRecordTag::PKMEAN_2,                 {"ps", true}  },
    { BamRecordTag::PKMID,                    {"pm", true}  },
    { BamRecordTag::PKMID_2,                  {"pi", true}  },
    { BamRecordTag::PRE_PULSE_FRAMES,         {"pd", true}  },
    { BamRecordTag::PULSE_CALL,               {"pc", true}  },
    { BamRecordTag::PULSE_CALL_WIDTH,         {"px", true}  },
    { BamRecordTag::PULSE_EXCLUSION,          {"pe", true}  },
    { BamRecordTag::PULSE_MERGE_QV,           {"pg", true}  },
    { BamRecordTag::PULSE_WIDTH,              {"pw", false} }, // 'pulse' in the name; but stored per-base, not per-pulse
    { BamRecordTag::REVERSE_IPD,              {"ri", false} },
    { BamRecordTag::REVERSE_PW,               {"rp", false} },
    { BamRecordTag::QUERY_END,                {"qe", false} },
    { BamRecordTag::QUERY_END_FRAME_NUMBER,   {"we", false} },
    { BamRecordTag::QUERY_START,              {"qs", false} },
    { BamRecordTag::QUERY_START_FRAME_NUMBER, {"ws", false} },
    { BamRecordTag::READ_ACCURACY,            {"rq", false} },
    { BamRecordTag::READ_GROUP,               {"RG", false} },
    { BamRecordTag::SCRAP_REGION_TYPE,        {"sc", false} },
    { BamRecordTag::SCRAP_ZMW_TYPE,           {"sz", false} },
    { BamRecordTag::SIGNAL_TO_NOISE,          {"sn", false} },
    { BamRecordTag::START_FRAME,              {"sf", true}  },
    { BamRecordTag::SUBSTITUTION_QV,          {"sq", false} },
    { BamRecordTag::SUBSTITUTION_TAG,         {"st", false} },
    { BamRecordTag::BASEMOD_LOCI,             {"MM", false} },
    { BamRecordTag::BASEMOD_QV,               {"ML", false} },

    // faux tags
    { BamRecordTag::SEQ,  {"  ",  false} },
    { BamRecordTag::QUAL, {"  ", false} }
};
// clang-format on

}  // namespace BAM
}  // namespace PacBio
