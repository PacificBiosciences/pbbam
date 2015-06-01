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

#ifndef PBIINDEX_H
#define PBIINDEX_H

#include "pbbam/Config.h"
#include "pbbam/PbiFile.h"
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal { class PbiIndexPrivate; }

// -----------------------------------------
//
// contiguous reads that satisfy a query will be returned as a ResultBlock
// this is to help minimize number of seeks (or even unneccesary checks)
//
// so an index query should expect 'ResultBlocks' as the return type
//
struct PBBAM_EXPORT IndexResultBlock
{
public:
    IndexResultBlock(int64_t o, size_t n);

public:
    int64_t offset_;
    size_t  numReads_;
};

typedef std::vector<IndexResultBlock> IndexResultBlocks;

class PBBAM_EXPORT PbiIndex
{
public:

public:
    /// \name Constructors & Related Methods
    /// \{

    PbiIndex(const std::string& pbiFilename);
    PbiIndex(const PbiIndex& other);
    PbiIndex(PbiIndex&& other);
    PbiIndex& operator=(const PbiIndex& other);
    PbiIndex& operator=(PbiIndex&& other);
    ~PbiIndex(void);

    /// \}

public:
    // PBI attributes
    bool HasBarcodeData(void) const;
    bool HasMappedData(void) const;
    bool HasReferenceData(void) const;
    bool HasSection(const PbiFile::Section section) const;

    PbiFile::Sections FileSections(void) const;
    uint32_t NumReads(void) const;
    PbiFile::VersionEnum Version(void) const;

public:
    IndexResultBlocks ReadGroupId(const std::string& value);

private:
    std::unique_ptr<internal::PbiIndexPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // PBIINDEX_H
