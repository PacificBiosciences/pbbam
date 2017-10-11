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
/// \file BamReader.h
/// \brief Defines the BamReader class.
//
// Author: Derek Barnett

#ifndef BAMREADER_H
#define BAMREADER_H

#include "pbbam/BamFile.h"
#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/GenomicInterval.h"

#include <htslib/sam.h>
#include <cstdint>
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

namespace internal { struct BamReaderPrivate; }

/// \brief The BamReader class provides basic read-access to a %BAM file.
///
/// The base-class implementation provides a sequential read-through of BAM
/// records. Derived classes may implement other access schemes (e.g. genomic
/// region, PBI-enabled record filtering).
///
class PBBAM_EXPORT BamReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Opens BAM file for reading.
    ///
    /// \param[in] fn %BAM filename
    /// \throws std::runtime_error if failed to open
    ///
    explicit BamReader(const std::string& fn);

    /// \brief Opens BAM file for reading.
    ///
    /// \param[in] bamFile BamFile object
    /// \throws std::runtime_error if failed to open
    ///
    explicit BamReader(const BamFile& bamFile);

    /// \brief Opens BAM file for reading.
    ///
    /// \param[in] bamFile BamFile object
    /// \throws std::runtime_error if failed to open
    ///
    explicit BamReader(BamFile&& bamFile);

    virtual ~BamReader();

    /// \}

public:
    /// \name BAM File Attributes
    /// \{

    /// \returns the underlying BamFile
    const BamFile& File() const;

    /// \returns %BAM filename
    std::string Filename() const;

    /// \returns BamHeader object from %BAM header contents
    const BamHeader& Header() const;

    /// \}

public:
    /// \name BAM File I/O
    /// \{

    /// \brief Fetches the "next" %BAM record.
    ///
    /// Default implementation will read records until EOF. Derived readers may
    /// use additional criteria to decide which record is "next" and when
    /// reading is done.
    ///
    /// \param[out] record  next BamRecord object. Should not be used if method
    ///                     returns false.
    ///
    /// \returns true if record was read successfully. Returns false if EOF (or
    ///          end of iterator in derived readers). False is not an error,
    ///          it indicates "end of data".
    ///
    /// \throws std::runtime_error if failed to read from file (e.g. possible
    ///         truncated or corrupted file).
    ///
    bool GetNext(BamRecord& record);

    /// \brief Seeks to virtual offset in %BAM.
    ///
    /// \note This is \b NOT a normal file offset, but the virtual offset used
    ///       in %BAM indexing.
    ///
    /// \throws std::runtime_error if failed to seek
    ///
    void VirtualSeek(int64_t virtualOffset);

    /// \returns current (virtual) file position.
    ///
    /// \note This is \b NOT a normal file offset, but the virtual offset used
    ///       in %BAM indexing.
    ///
    int64_t VirtualTell() const;

    /// \}

protected:
    /// \name BAM File I/O
    /// \{

    /// \brief Helper method for access to underlying BGZF stream pointer.
    ///
    /// Useful for derived readers' contact points with htslib methods.
    ///
    /// \returns BGZF stream pointer
    ///
    BGZF* Bgzf() const;

    /// \brief Performs the actual raw read of the next record from the BAM
    ///        file.
    ///
    /// Default implementation will read records, sequentially, until EOF.
    /// Derived readers may use additional criteria to decide which record is
    ///  "next" and when reading is done.
    ///
    /// Return value should be equivalent to htslib's bam_read1():
    ///     >= 0 : normal
    ///       -1 : EOF (not an error)
    ///     < -1 : error
    ///
    /// \param[in]  bgzf BGZF stream pointer
    /// \param[out] b    %BAM record pointer
    /// \returns integer status code, see description
    ///
    virtual int ReadRawData(BGZF* bgzf, bam1_t* b);

    /// \}

private:
    std::unique_ptr<internal::BamReaderPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // BAMREADER_H
