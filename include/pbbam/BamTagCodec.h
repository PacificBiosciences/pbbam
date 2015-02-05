#ifndef BAMTAGCODEC_H
#define BAMTAGCODEC_H

#include "pbbam/Config.h"
#include "pbbam/TagCollection.h"
#include <vector>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT BamTagCodec
{

// high-level, operate on a full collection
public:
    static TagCollection Decode(const std::vector<uint8_t>& data);
    static std::vector<uint8_t> Encode(const PacBio::BAM::TagCollection& tags);

// per-tag methods
public:

    // returns the SAM/BAM single char code for tag type
    static uint8_t TagTypeCode(const PacBio::BAM::Tag& tag);

    // returns the tag value's raw data in bytes (does *NOT* encode name & type)
    static std::vector<uint8_t> ToRawData(const PacBio::BAM::Tag& tag);

    // TODO: make this hidden a bit more, maybe this whole class in fact
    // rawData should be the result of sam.h:bam_aux_get(...)
    static PacBio::BAM::Tag FromRawData(uint8_t* rawData);
};

} // namespace BAM
} // namespace PacBio

#endif // BAMTAGCODEC_H
