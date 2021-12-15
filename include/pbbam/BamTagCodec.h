#ifndef PBBAM_BAMTAGCODEC_H
#define PBBAM_BAMTAGCODEC_H

#include <pbbam/Config.h>

#include <pbbam/TagCollection.h>

#include <vector>

#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief The BamTagCodec class provides binary encoding/decoding of %BAM tag
///        data.
///
/// \note BamTagCodec is mostly an implementation and/or testing detail, and may
///       be removed from the public API.
///
class PBBAM_EXPORT BamTagCodec
{
public:
    /// \name Tag Collection Methods
    /// \{

    /// \brief Creates a TagCollection from raw BAM data.
    ///
    /// \param[in] data     BAM-formatted (binary) tag data
    /// \returns TagCollection containing tag data
    ///
    static TagCollection Decode(const std::vector<uint8_t>& data);

    /// \brief Creates binary BAM data from a TagCollection.
    ///
    /// \param[in] tags     TagCollection containing tag data
    /// \returns vector of bytes (encoded BAM data)
    ///
    static std::vector<uint8_t> Encode(const TagCollection& tags);

    /// \}

public:
    /// \name Per-Tag Methods
    /// \{

    /// \brief Determines the SAM/BAM tag code for a Tag.
    ///
    /// \param[in] tag                  Tag object to check
    /// \param[in] additionalModifier   optional extra modifier (allows explicit
    ///                                 modification of an otherwise const Tag)
    ///
    /// \returns the SAM/BAM single char code for tag type
    ///
    static uint8_t TagTypeCode(const Tag& tag,
                               const TagModifier& additionalModifier = TagModifier::NONE);

    /// \brief Encodes a single Tag's contents in %BAM binary
    ///
    /// \note This method does \b NOT encode the tag name & tag type. It does
    ///       include the element type for array-type tags.
    ///
    /// \param[in] tag                  Tag object containing data to encode
    /// \param[in] additionalModifier   optional extra modifier (allows explicit
    ///                                 modification of an otherwise const Tag)
    ///
    /// \returns vector of bytes (encoded BAM data)
    ///
    static std::vector<uint8_t> ToRawData(
        const Tag& tag, const TagModifier& additionalModifier = TagModifier::NONE);

    /// \brief Creates a Tag object from binary BAM data.
    ///
    /// \param[in] rawData      raw BAM bytes (assumed to be the result of
    ///                         htslib's bam_aux_get())
    ///
    /// \returns resulting Tag object
    ///
    static Tag FromRawData(uint8_t* rawData);

    /// \}
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BAMTAGCODEC_H
