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
/// \file SamWriter.h
/// \brief Defines the SamWriter class.
//
// Author: Derek Barnett

#ifndef SAMWRITER_H
#define SAMWRITER_H

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/IRecordWriter.h"
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

namespace internal { class SamWriterPrivate; }

/// \brief The SamWriter class provides a writing interface for creating
///        new SAM files.
///
/// \note The underlying buffered data may not be flushed to the file until the
///       destructor is called. Trying to access the file (reading, stat-ing,
///       indexing, etc.) before the SamWriter is destroyed yields undefined
///       behavior. Enclose the SamWriter in some form of local scope (curly
///       braces, a separate function, etc.) to ensure that its destructor is
///       called before proceeding to read-based operations.
///
/// \code{.cpp}
///  {
///     SamWriter w(...);
///     // write data
///  }
///  // now safe to access the new file
/// \endcode
///
///
class SamWriter : public IRecordWriter
{
public:
    /// \brief Opens a SAM file for writing & writes the header information.
    ///
    /// \note Set \p filename to "-" for stdout.
    ///
    /// \param[in] filename     path to output SAM file
    /// \param[in] header       BamHeader object
    ///
    /// \throws std::runtime_error if there was a problem opening the file for
    ///         writing or if an error occurred while writing the header
    ///
    SamWriter(const std::string& filename, const BamHeader& header);

    /// Fully flushes all buffered data & closes file.
    ///
    ~SamWriter() override;

    /// Copy and Move constructors are disabled
    SamWriter(const SamWriter&) = delete;
    SamWriter(SamWriter&&) = delete;
    SamWriter& operator=(const SamWriter&) = delete;
    SamWriter& operator=(SamWriter&&) = delete;

public:
    
    /// \brief Try to flush any buffered data to file.
    ///
    /// \note The underlying implementation may not necessarily flush buffered
    ///       data immediately, especially in a multithreaded writer situation.
    ///       Let the SamWriter go out of scope to fully ensure flushing.
    ///
    /// \throws std::runtime_error if flush fails
    ///
    void TryFlush() override;

    /// \brief Write a record to the output SAM file.
    ///
    /// \param[in] record BamRecord object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    void Write(const BamRecord& record) override;

    /// \brief Write a record to the output SAM file.
    ///
    /// \param[in] recordImpl BamRecordImpl object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    void Write(const BamRecordImpl& recordImpl) override;

private:
    std::unique_ptr<internal::SamWriterPrivate> d_;
};

} // namesapce BAM
} // namespace PacBio

#endif // SAMWRITER_H
