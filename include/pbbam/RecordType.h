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
inline bool IsCcsOrTranscript(const RecordType type)
{
    return type == RecordType::CCS || type == RecordType::TRANSCRIPT;
}

}  // namespace BAM
}  // namespace PacBio

#endif  // RECORDTYPE_H
