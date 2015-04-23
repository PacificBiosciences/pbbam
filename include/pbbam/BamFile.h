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

namespace internal { class BamFilePrivate; }

class PBBAM_EXPORT BamFile
{
public:

    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a BamFile object on the provided \p filename & loads header information.
    ///
    /// \param[in] filename BAM filename
    /// \throws std::exception on failure
    BamFile(const std::string& filename);

    BamFile(const BamFile& other);
    BamFile(BamFile&& other);
    BamFile& operator=(const BamFile& other);
    BamFile& operator=(BamFile&& other);
    ~BamFile(void);

    /// \}

public:

    /// \name Filename Methods
    /// \{

    /// \returns BAM filename
    std::string Filename(void) const;

    /// \returns filename of PacBio index file (".pbi")
    /// \note No guarantee is made on the existence of this file.
    ///       This method simply returns the expected filename.
    std::string PacBioIndexFilename(void) const;

    /// \returns filename of standard index file (".bai")
    /// \note No guarantee is made on the existence of this file.
    ///       This method simply returns the expected filename.
    std::string StandardIndexFilename(void) const;

    /// \}

    /// \name Header Metadata Methods
    /// \{

    /// \returns true if header metadata has this reference name
    bool HasReference(const std::string& name) const;

    /// \returns BamHeader containing the file's metadata
    BamHeader Header(void) const;

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

private:
    PBBAM_SHARED_PTR<internal::BamFilePrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // BAMFILE_H
