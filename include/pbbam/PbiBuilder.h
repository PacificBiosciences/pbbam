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
// Author: Derek Barnett

#ifndef PBIBUILDER_H
#define PBIBUILDER_H

#include "pbbam/Config.h"
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

class BamRecord;
class PbiRawData;

namespace internal { class PbiBuilderPrivate; }

/// This class may be used to construct PBI index data while a BAM file is being
/// written, rather than waiting to process it at the end.
///
class PBBAM_EXPORT PbiBuilder
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// Initialize builder to write data to \p pbiFilename.
    ///
    /// \throws std::runtime_error if PBI file cannot be opened for writing
    ///
    PbiBuilder(const std::string& pbiFilename);

    /// Initialize builder to write data to \p pbiFilename. Reference data-tracking
    /// structures will be initialized to expect \p numReferenceSequences. (This is
    /// useful so that we can mark any references that lack observed data appropriately).
    ///
    /// \throws std::runtime_error if PBI file cannot be opened for writing
    ///
    PbiBuilder(const std::string& pbiFilename, const size_t numReferenceSequences);

    /// On destruction, data summaries are calculated, raw data is written to file, and
    /// file handle closed.
    ///
    ~PbiBuilder(void);

    /// \}

public:
    /// \name Index Building
    /// \{

    /// Adds \p record's data to underlying raw data structure. \p vOffset is the BGZF
    /// virtual offset into the BAM file where the record begins.
    ///
    /// \sa BamWriter::Write(const BamRecord& record, int64_t* vOffset) for the easiest
    ///     way to retrieve this information while generating a BAM file. See example below:
    ///
    /// \code{.cpp}
    ///  BamWriter writer(...);
    ///  PbiBuilder pbiBuilder(...);
    ///  int64_t vOffset;
    ///  while (...) {
    ///     BamRecord record;
    ///     // ... generate record data ...
    ///     writer.Write(record, &vOffset);
    ///     pbiBuilder.AddRecord(record, &vOffset);
    ///  }
    /// \endcode
    ///
    void AddRecord(const BamRecord& record, const int64_t vOffset);

    /// \returns const reference to current raw index data. Mostly only used for testing;
    ///          shouldn't be needed by most client code.
    ///
    const PbiRawData& Index(void) const;

    /// \}

private:
    std::unique_ptr<internal::PbiBuilderPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // PBIBUILDER_H
