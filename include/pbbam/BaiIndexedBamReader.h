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
/// \file BaiIndexedBamReader.h
/// \brief Defines the BaiIndexedBamReader class.
//
// Author: Derek Barnett

#ifndef BAIINDEXEDBAMREADER_H
#define BAIINDEXEDBAMREADER_H

#include "pbbam/BamReader.h"
#include "pbbam/BamFile.h"
#include "pbbam/GenomicInterval.h"

namespace PacBio {
namespace BAM {

namespace internal { struct BaiIndexedBamReaderPrivate; }

/// \brief The BaiIndexedBamReader class provides read-only iteration over %BAM
///        records, bounded by a particular genomic interval.
///
/// The SAM/BAM standard index (*.bai) is used to allow random-access operations.
///
class PBBAM_EXPORT BaiIndexedBamReader : public BamReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Constructs %BAM reader, bounded by a genomic interval.
    ///
    /// All reads that overlap the interval will be available.
    ///
    /// \param[in] interval iteration will be bounded by this GenomicInterval.
    /// \param[in] filename input %BAM filename
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    BaiIndexedBamReader(const GenomicInterval& interval,
                        const std::string& filename);

    /// \brief Constructs BAM reader, bounded by a genomic interval.
    ///
    /// All reads that overlap the interval will be available.
    ///
    /// \param[in] interval iteration will be bounded by this GenomicInterval.
    /// \param[in] bamFile input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    BaiIndexedBamReader(const GenomicInterval& interval, const BamFile& bamFile);

    /// \brief Constructs %BAM reader, bounded by a genomic interval.
    ///
    /// All reads that overlap the interval will be available.
    ///
    /// \param[in] interval iteration will be bounded by this GenomicInterval.
    /// \param[in] bamFile input BamFile object
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    BaiIndexedBamReader(const GenomicInterval& interval, BamFile&& bamFile);

    /// \}

public:
    /// \name Random-Access
    /// \{

    /// \returns the current GenomicInterval in use by this reader
    const GenomicInterval& Interval() const;

    /// \brief Sets a new genomic interval on the reader.
    ///
    /// \param[in] interval
    /// \returns reference to this reader
    ///
    BaiIndexedBamReader& Interval(const GenomicInterval& interval);

    /// \}

protected:
    int ReadRawData(BGZF* bgzf, bam1_t* b) override;

private:
    std::unique_ptr<internal::BaiIndexedBamReaderPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // BAIINDEXEDBAMREADER_H
