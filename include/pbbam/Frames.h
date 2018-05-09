// File Description
/// \file Frames.h
/// \brief Defines the Frames class.
//
// Author: Derek Barnett

#ifndef FRAMES_H
#define FRAMES_H

#include <cstddef>
#include <cstdint>
#include <vector>
#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief The Frames class represents pulse frame data.
///
/// Frame data may be stored in either their raw, 16-bit values or
/// using a lossy, 8-bit compression scheme.
///
/// This class is used to store the data and convert between the 2 storage types.
///
class PBBAM_EXPORT Frames
{
public:
    /// \name Conversion Methods
    /// \{

    /// \brief Constructs a Frames object from encoded (lossy, 8-bit) data.
    ///
    /// \note This method should probably not be needed often by client code
    ///       working with frame data. It exists primarily for (internal)
    ///       parsing & interpretation of the %BAM file contents. The method is
    ///       available, though, should the conversion operation be needed.
    ///
    /// \param[in] codedData    encoded data
    /// \returns Frames object
    ///
    static Frames Decode(const std::vector<uint8_t>& codedData);

    /// \brief Creates encoded, compressed frame data from raw input data.
    ///
    /// \param[in] frames   raw frame data
    /// \returns lossy, 8-bit encoded frame data
    ///
    static std::vector<uint8_t> Encode(const std::vector<uint16_t>& frames);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    Frames(std::vector<uint16_t> frames);

    Frames();
    Frames(const Frames&) = default;
    Frames(Frames&&) = default;
    Frames& operator=(const Frames&) = default;
    Frames& operator=(Frames&&) = default;
    ~Frames() = default;

    /// \}

public:
    /// \name Access Data
    /// \{

    /// \returns Frame data in expanded (not encoded) form
    std::vector<uint16_t>& DataRaw();
    const std::vector<uint16_t>& Data() const;

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// \returns Frame data in (lossy, 8-bit) encoded form.
    std::vector<uint8_t> Encode() const;

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    bool operator==(const Frames& other) const;
    bool operator!=(const Frames& other) const;

    /// \}

public:
    /// \name STL Compatbility
    /// \{

    /// \returns A const_iterator to the beginning of the sequence.
    std::vector<uint16_t>::const_iterator cbegin() const;

    /// \returns A const_iterator to the element past the end of the sequence.
    std::vector<uint16_t>::const_iterator cend() const;

    /// \returns A const_iterator to the beginning of the sequence.
    std::vector<uint16_t>::const_iterator begin() const;

    /// \returns A const_iterator to the element past the end of the sequence.
    std::vector<uint16_t>::const_iterator end() const;

    /// \returns An iterator to the beginning of the sequence.
    std::vector<uint16_t>::iterator begin();

    /// \returns An iterator to the element past the end of the sequence.
    std::vector<uint16_t>::iterator end();

    /// \returns The number of frame data points.
    size_t size() const;

    /// \returns True if the container is empty, false otherwise.
    bool empty() const;

    /// \}

public:
    /// \name Access Data
    /// \{

    /// Sets this record's data.
    ///
    /// \param[in] frames data in expanded (not encoded) form
    /// \returns reference to this object
    ///
    Frames& Data(std::vector<uint16_t> frames);

    /// \}

private:
    std::vector<uint16_t> data_;
};

}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/Frames.inl"

#endif  // FRAMES_H
