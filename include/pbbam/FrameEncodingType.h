// File Description
/// \file FrameEncodingType.h
/// \brief Defines the FrameEncodingType enum.
//
// Author: Derek Barnett

#ifndef FRAMEENCODINGTYPE_H
#define FRAMEENCODINGTYPE_H

namespace PacBio {
namespace BAM {

/// \brief This enum defines the possible encoding modes used in Frames data
/// (e.g. BamRecord::IPD or BamRecord::PulseWidth).
///
/// The LOSSY mode is the default in production output; LOSSLESS mode
/// being used primarily for internal applications.
///
/// \sa https://github.com/PacificBiosciences/PacBioFileFormats/blob/3.0/BAM.rst
///     for more information on pulse frame encoding schemes.
///
enum class FrameEncodingType
{
    LOSSY,    ///< 8-bit compression (using CodecV1) of frame data
    LOSSLESS  ///< 16-bit native frame data
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FRAMEENCODINGTYPE_H
