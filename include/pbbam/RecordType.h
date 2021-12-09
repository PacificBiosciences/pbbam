#ifndef PBBAM_RECORDTYPE_H
#define PBBAM_RECORDTYPE_H

#include <pbbam/Config.h>

#include <string>

namespace PacBio {
namespace BAM {

/// \brief This enum defines the possible PacBio BAM record types.
///
/// \sa ReadGroupInfo::ReadType
///
enum class RecordType
{
    ZMW,         ///< Polymerase read
    HQREGION,    ///< High-quality region
    SUBREAD,     ///< Subread
    CCS,         ///< Circular consensus sequence
    SCRAP,       ///< Additional sequence (barcodes, adapters, etc.)
    UNKNOWN,     ///< Unknown read type
    TRANSCRIPT,  ///< Transcript

    POLYMERASE = ZMW  ///< \deprecated as of PacBio BAM spec v 3.0.4 (use RecordType::ZMW instead)
};

///
/// \brief IsCcsOrTranscript
///
/// CCS & Transcript type records handle queryStart/End in the same way. This
/// status is checked in several places, so this is a convenient helper.
///
/// \param[in] type
///
bool IsCcsOrTranscript(RecordType type);

///
/// \brief Returns string representation of RecordType
///
/// \param type
/// \return std::string
/// \throws std::runtime_error if type is unrecognized
///
std::string ToString(RecordType type);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_RECORDTYPE_H
