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
/// \file PbiIndexedBamReader.h
/// \brief Defines the PbiIndexedBamReader class.
//
// Author: Derek Barnett

#ifndef PBIINDEXEDBAMREADER_H
#define PBIINDEXEDBAMREADER_H

#include "pbbam/BamFile.h"
#include "pbbam/BamReader.h"
#include "pbbam/PbiBasicTypes.h"
#include "pbbam/PbiFilter.h"
#include "pbbam/PbiIndex.h"
#include <string>

namespace PacBio {
namespace BAM {

namespace internal { struct PbiIndexedBamReaderPrivate; }

/// \brief The PbiIndexedBamReader class provides read-only iteration over %BAM
///        records, limited to some filtering criteria.
///
/// The PacBio BAM index (*.pbi) is used to allow random-access operations.
///
class PBBAM_EXPORT PbiIndexedBamReader : public BamReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Constructs %BAM reader, with an initial filter.
    ///
    /// All reads that satisfy the filter will be available.
    ///
    /// \param[in] filter       PbiFilter or compatible object
    /// \param[in] bamFilename  input %BAM filename
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(const PbiFilter& filter, const std::string& bamFilename);

    /// \brief Constructs %BAM reader, with an initial filter.
    ///
    /// All reads that satisfy the filter will be available.
    ///
    /// \param[in] filter       PbiFilter or compatible object
    /// \param[in] bamFile      input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(const PbiFilter& filter, const BamFile& bamFile);

    /// \brief Constructs %BAM reader, with an initial filter.
    ///
    /// All reads that satisfy the filter will be available.
    ///
    /// \param[in] filter       PbiFilter or compatible object
    /// \param[in] bamFile      input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(const PbiFilter& filter, BamFile&& bamFile);

    /// \brief Constructs %BAM reader, with no initial filter.
    ///
    /// Useful for delaying either specifying the filtering criteria or
    /// performing the PBI lookups.
    ///
    /// \param[in] bamFilename  input %BAM filename
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(const std::string& bamFilename);

    /// \brief Constructs %BAM reader, with no initial filter.
    ///
    /// Useful for delaying either specifying the filtering criteria or
    /// performing the PBI lookups.
    ///
    /// \param[in] bamFile      input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(const BamFile& bamFile);

    /// \brief Constructs %BAM reader, with no initial filter.
    ///
    /// Useful for delaying either specifying the filtering criteria or
    /// performing the PBI lookups.
    ///
    /// \param[in] bamFile      input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.pbi) cannot be
    ///         read
    ///
    PbiIndexedBamReader(BamFile&& bamFile);

    ~PbiIndexedBamReader() override;

    /// \}

public:
    /// \name Filtering & Index Data
    /// \{

    /// \returns the current filter active on this reader
    const PbiFilter& Filter() const;

    uint32_t NumReads() const;

//    /// \returns the reader's underlying index data
//    const PbiIndex& Index() const;

public:
    /// \brief Sets a new filter on the reader.
    ///
    /// \param[in] filter
    /// \returns reference to this reader
    ///
    PbiIndexedBamReader& Filter(const PbiFilter& filter);

    /// \}

protected:
    int ReadRawData(BGZF* bgzf, bam1_t* b) override;

private:
    std::unique_ptr<internal::PbiIndexedBamReaderPrivate> d_;
};

} // namespace internal
} // namespace BAM

#endif // PBIINDEXEDBAMREADER_H
