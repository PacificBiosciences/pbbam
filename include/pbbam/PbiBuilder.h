// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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
/// \file PbiBuilder.h
/// \brief Defines the PbiBuilder class.
//
// Author: Derek Barnett

#ifndef PBIBUILDER_H
#define PBIBUILDER_H

#include "pbbam/Config.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

class BamRecord;
class PbiRawData;

namespace internal { class PbiBuilderPrivate; }

/// \brief The PbiBuilder class construct PBI index data from %BAM record data.
///
/// Records are added one-by-one. This allows for either whole-file indexing of
/// existing %BAM files or for indexing "on-the-fly" alongside a %BAM file as it
/// is generated.
///
/// For simple PBI creation from existing %BAM files, see PbiFile::CreateFrom.
/// This is the recommended approach, unless finer control or additional
/// processing is needed.
///
class PBBAM_EXPORT PbiBuilder
{
public:
    /// \brief This enum allows you to control the compression level of the
    ///        output PBI file.
    ///
    /// Values are equivalent to zlib compression levels. See its documentation
    /// for more details: http://www.zlib.net/manual.html
    ///
    enum CompressionLevel
    {
        CompressionLevel_0 = 0
      , CompressionLevel_1 = 1
      , CompressionLevel_2 = 2
      , CompressionLevel_3 = 3
      , CompressionLevel_4 = 4
      , CompressionLevel_5 = 5
      , CompressionLevel_6 = 6
      , CompressionLevel_7 = 7
      , CompressionLevel_8 = 8
      , CompressionLevel_9 = 9

      , DefaultCompression = -1
      , NoCompression      = CompressionLevel_0
      , FastCompression    = CompressionLevel_1
      , BestCompression    = CompressionLevel_9
    };

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Initializes builder to write data to \p pbiFilename.
    ///
    /// \param[in] pbiFilename      output filename
    /// \param[in] compressionLevel zlib compression level
    /// \param[in] numThreads       number of threads for compression. If set to
    ///                             0, PbiBuilder will attempt to determine a
    ///                             reasonable estimate. If set to 1, this will
    ///                             force single-threaded execution. No checks
    ///                             are made against an upper limit.
    ///
    /// \throws std::runtime_error if PBI file cannot be opened for writing
    ///
    PbiBuilder(const std::string& pbiFilename,
               const PbiBuilder::CompressionLevel compressionLevel = PbiBuilder::DefaultCompression,
               const size_t numThreads = 4);

    /// \brief Initializes builder to write data to \p pbiFilename.
    ///
    /// Reference data-tracking structures will be initialized to expect
    /// \p numReferenceSequences. (This is useful so that we can mark any
    /// references that lack observed data appropriately).
    ///
    /// \param[in] pbiFilename              output filename
    /// \param[in] numReferenceSequences    number of possible reference
    ///                                     sequences, e.g. BamHeader::NumSequences
    /// \param[in] compressionLevel zlib compression level
    /// \param[in] numThreads       number of threads for compression. If set to
    ///                             0, PbiBuilder will attempt to determine a
    ///                             reasonable estimate. If set to 1, this will
    ///                             force single-threaded execution. No checks
    ///                             are made against an upper limit.
    ///
    /// \throws std::runtime_error if PBI file cannot be opened for writing
    ///
    PbiBuilder(const std::string& pbiFilename,
               const size_t numReferenceSequences,
               const PbiBuilder::CompressionLevel compressionLevel = PbiBuilder::DefaultCompression,
               const size_t numThreads = 4);

    /// \brief Initializes builder to write data to \p pbiFilename.
    ///
    /// Reference data-tracking structures will be initialized to expect
    /// \p numReferenceSequences, but only if \p isCoordinateSorted is true.
    ///
    /// \param[in] pbiFilename              output filename
    /// \param[in] numReferenceSequences    number of possible reference
    ///                                     sequences, e.g. BamHeader::NumSequences
    /// \param[in] isCoordinateSorted       if false, disables reference
    ///                                     sequence tracking
    ///                                     (BamHeader::SortOrder != "coordinate")
    /// \param[in] compressionLevel zlib compression level
    /// \param[in] numThreads       number of threads for compression. If set to
    ///                             0, PbiBuilder will attempt to determine a
    ///                             reasonable estimate. If set to 1, this will
    ///                             force single-threaded execution. No checks
    ///                             are made against an upper limit.
    ///
    /// \throws std::runtime_error if PBI file cannot be opened for writing
    ///
    PbiBuilder(const std::string& pbiFilename,
               const size_t numReferenceSequences,
               const bool isCoordinateSorted,
               const PbiBuilder::CompressionLevel compressionLevel = PbiBuilder::DefaultCompression,
               const size_t numThreads = 4);

    /// \brief Destroys builder, writing its data out to PBI file.
    ///
    ///
    /// \note Exceptions are swallowed. Use Close() if you want to catch them.
    ///
    ~PbiBuilder() noexcept;

    /// \}

public:
    /// \name Index Building
    /// \{

    /// \brief Adds \p record's data to underlying raw data structure.
    ///
    /// \note \p vOffset is a BGZF \b virtual offset into the %BAM file. To get
    ///          this value, you should use one of the following: \n
    ///        - while reading existing %BAM: BamReader::VirtualTell \n
    ///        - while writing new %BAM:      BamWriter::Write(const BamRecord& record, int64_t* vOffset) \n
    ///
    ///
    /// To build a PBI index while generating a %BAM file:
    /// \include code/PbiBuilder_WithWriter.txt
    ///
    /// To build a PBI index from an existing %BAM file:
    /// \include code/PbiBuilder_WithReader.txt
    ///
    /// \param[in] record   input BamRecord to pull index data from
    /// \param[in] vOffset  \b virtual offset into %BAM file where record begins
    ///
    void AddRecord(const BamRecord& record, const int64_t vOffset);

    /// \brief Writes data out to PBI file & closes builder.
    ///
    /// \note Any exceptions are thrown to caller. If you don't care about
    ///       catching exceptions with file I/O, just let the builder go out of
    ///       scope and data will be written, but exceptions swallowed (to avoid
    ///       throwing from destructor).
    ///
    void Close();

    /// \}

private:
    std::unique_ptr<internal::PbiBuilderPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // PBIBUILDER_H
