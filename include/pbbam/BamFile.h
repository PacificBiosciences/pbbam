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

// Author: Derek Barnett

#ifndef BAMFILE_H
#define BAMFILE_H

#include "pbbam/Config.h"
#include "pbbam/BamHeader.h"
#include <string>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT BamFile
{
public:

    /// This enum describes the errors that may be returned by the Error() function.
    enum FileError
    {
         NoError         ///< No error occurred.
       , OpenError       ///< An error occurred when opening the file.
       , ReadHeaderError ///< An error occurred when reading the header data.
    };

public:

    /// \name Constructors & Related Methods
    /// \{

    /// Creates a BamFile object with no associated file.
    BamFile(void);

    /// \brief Creates a BamFile object on the provided \p filename.
    ///
    /// \sa Open
    /// \param[in] filename BAM filename
    BamFile(const std::string& filename);

    /// Creates a copy of the BamFile object, with the same filename, header data, and error status
    BamFile(const BamFile& other);

    ~BamFile(void);

    /// \}

public:

    /// \name Filename Methods
    /// \{

    /// \returns BAM filename
    std::string Filename(void) const;

    /// \returns filename of PacBio index file (".pbi")
    /// \note No guarantee is made on the existence of this file - this method simply returns the
    ///       expected filename.
    std::string PacBioIndexFilename(void) const;

    /// \returns filename of standard index file (".bai")
    /// \note No guarantee is made on the existence of this file - this method simply returns the
    ///       expected filename.
    std::string StandardIndexFilename(void) const;

    /// \}

    /// \name Header Metadata Methods
    /// \{

    /// \returns BamHeader containing the file's metadata
    BamHeader::SharedPtr Header(void) const;

    /// \returns true if BAM file is a PacBio BAM file (i.e. has non-empty version associated with header "pb" tag)
    bool IsPacBioBAM(void) const;

    /// \returns ID for reference \p name (can be used for e.g. GenomicIntervalQuery), -1 if not found
    int ReferenceId(const std::string& name) const;

    /// \return name of reference matching \p id, empty string if not found
    std::string ReferenceName(const int id) const;

    /// \returns length of requested reference \p name. 0 if not found
    uint32_t ReferenceLength(const std::string& name) const;

    /// \returns length of requested reference \p id. 0 if not found
    uint32_t ReferenceLength(const int id) const;

    /// \}

public:
    /// \name Open/Close Methods
    /// \{

    /// \returns true if BamFile has been opened on a file
    bool IsOpen(void) const;

    /// \}

public:
    /// \name Error Handling
    /// \{

    /// \returns file error status
    BamFile::FileError Error(void) const;

    /// \returns true if BamFile::Error() == NoError
    ///
    /// \code
    ///     BamFile file("foo.bam");
    ///     if (!file) {
    ///         // handle error
    ///         return;
    ///     }
    ///     // ok to work with file and/or header info
    /// \endcode
    operator bool(void) const;

    /// \}

public:
    /// \name Open/Close Methods
    /// \{

    /// Resets BAM file metadata associated with this object.
    void Close(void);

    /// Attempts to open the file and load the header metadata. The error status
    /// will be set if either operation fails.
    void Open(const std::string& filename);

    /// \}

private:
    std::string filename_;
    BamFile::FileError error_;
    BamHeader::SharedPtr header_;
};

inline BamFile::operator bool(void) const
{ return error_ == BamFile::NoError; }

inline BamFile::FileError BamFile::Error(void) const
{ return error_; }

inline std::string BamFile::Filename(void) const
{ return filename_; }

inline BamHeader::SharedPtr BamFile::Header(void) const
{ return header_; }

inline bool BamFile::IsPacBioBAM(void) const
{ return (header_ ? !header_->PacBioBamVersion().empty() : false); }

inline std::string BamFile::StandardIndexFilename(void) const
{ return filename_ + ".bai"; }

inline std::string BamFile::PacBioIndexFilename(void) const
{ return filename_ + ".pbi"; }

inline int BamFile::ReferenceId(const std::string& name) const
{ return (header_ ? header_->SequenceId(name) : -1); }

inline uint32_t BamFile::ReferenceLength(const std::string& name) const
{ return ReferenceLength(ReferenceId(name)); }

inline uint32_t BamFile::ReferenceLength(const int id) const
{ return (header_ ? std::stoul(header_->SequenceLength(id)) : 0); }

inline std::string BamFile::ReferenceName(const int id) const
{ return (header_ ? header_->SequenceName(id) : ""); }

} // namespace BAM
} // namespace PacBio

#endif // BAMFILE_H
