// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#ifndef BAMWRITER_H
#define BAMWRITER_H


#include "pbbam/Config.h"
#include "pbbam/SamHeader.h"
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

class BamFile;
class BamRecord;

class PBBAM_EXPORT BamWriter
{
public:
    /// This enum allows you to control the compression level of the output BAM file.
    ///
    /// Values are equivalent to zlib compression levels. See its documentation for more details:
    /// http://www.zlib.net/manual.html
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

    /// This enum describes the errors that may be returned by the Error() function.
    enum WriteError
    {
        NoError = 0      ///< No error occurred.
      , OpenFileError    ///< An error occurred while opening the file.
      , NullHeaderError  ///< Header data was invalid.
      , WriteHeaderError ///< An error occurred while writing header data.
      , WriteRecordError ///< An error occurred while writing a record.
    };

public:

    /// \name Constructors & Related Methods
    /// \{

    /// Opens a BAM file for writing & writes the header information.
    ///
    /// The error status will be set if either operation fails.
    ///
    /// \note Set \p filename to "-" for stdout.
    ///
    /// \param[in] filename path to output BAM file
    /// \param[in] header SamHeader object
    /// \param[in] compressionLevel zlib compression level
    BamWriter(const std::string& filename,
              const SamHeader& header,
              const BamWriter::CompressionLevel compressionLevel = BamWriter::DefaultCompression);

    ~BamWriter(void);

    /// \}

public:

    /// \name Data Writing & Resource Management

    /// Closes BAM file - flushing any buffered data & releasing the file handle.
    ///
    /// The destructor handles cleanup and should suffice in most cases. However, this method
    /// allows an earlier cleanup, if  necessary.
    void Close(void);

    /// Write any buffered data to file.
    ///
    /// \returns true if successful
    bool Flush(void);

    /// Write a record to the output BAM file.
    ///
    /// \param[in] record BamRecord object
    ///
    /// \returns succcess/failure
    bool Write(const BamRecord& record);

    /// \}

    /// \name Error Handling
    /// \{

    /// \returns error status code
    BamWriter::WriteError Error(void) const;

    /// \returns true if BamWriter::Error() == NoError
    ///
    /// \code
    ///     BamWriter writer(...);
    ///     if (!writer) {
    ///         // handle error
    ///         return;
    ///     }
    ///     // ok to work with writer
    /// \endcode
    operator bool(void) const;

    /// \}

private:
    bool Open(const std::string& filename,
              const std::shared_ptr<bam_hdr_t> rawHeader,
              const BamWriter::CompressionLevel compressionLevel = BamWriter::DefaultCompression);
    bool Write(const std::shared_ptr<bam1_t>& rawRecord);

private:
    std::shared_ptr<samFile>   file_;
    std::shared_ptr<bam_hdr_t> header_;
    std::string filename_;
    BamWriter::WriteError error_;
};

} // namespace BAM
} // namespace PacBio

#endif // BAMWRITER_H
