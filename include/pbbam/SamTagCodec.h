// File Description
/// \file SamTagCodec.h
/// \brief Defines the SamTagCodec class.
//
// Author: Derek Barnett

#ifndef SAMTAGCODEC_H
#define SAMTAGCODEC_H

#include <string>
#include "pbbam/Config.h"
#include "pbbam/TagCollection.h"

namespace PacBio {
namespace BAM {

/// \brief The SamTagCodec class provides text-based encoding/decoding of %BAM
///        tag data.
///
/// \note SamTagCodec is mostly an implementation and/or testing detail, and may
///       be removed from the public API.
///
class PBBAM_EXPORT SamTagCodec
{
public:
    /// \name Tag Collection Methods
    /// \{

    /// \brief Creates a TagCollection from SAM-formatted tag data.
    ///
    /// \param[in] tagString    SAM-formmated string
    /// \returns resulting tag collection
    ///
    static TagCollection Decode(const std::string& tagString);

    /// \brief Creates SAM-formatted string from a Tag.
    ///
    /// \param[in] name 2-character tag name
    /// \param[in] tag  Tag instance containing data
    ///
    /// \return SAM-formatted string
    ///
    static std::string Encode(const std::string& name, const PacBio::BAM::Tag& tag);

    /// \brief Creates SAM-formatted string from a TagCollection.
    ///
    /// \param[in] tags     TagCollection containing tag data
    /// \returns SAM-formatted string
    ///
    static std::string Encode(const PacBio::BAM::TagCollection& tags);
};

///
/// \brief creates a tag per the SAM/BAM text format
///
/// \param tag    tag name
/// \param value  tag value
///
/// \return formatted tag string
///
std::string MakeSamTag(std::string tag, std::string value);

}  // namespace BAM
}  // namespace PacBio

#endif  // SAMTAGCODEC_H
