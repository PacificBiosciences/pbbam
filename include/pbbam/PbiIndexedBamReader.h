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

class PBBAM_EXPORT PbiIndexedBamReader : public BamReader
{
public:
    PbiIndexedBamReader(const PbiFilter& filter, const std::string& bamFilename);
    PbiIndexedBamReader(const PbiFilter& filter, const BamFile& bamFile);
    PbiIndexedBamReader(const PbiFilter& filter, BamFile&& bamFile);

public:
    const PbiFilter& Filter(void) const;
    PbiIndexedBamReader& Filter(const PbiFilter& filter);

protected:
    int ReadRawData(BGZF* bgzf, bam1_t* b);

private:
    std::unique_ptr<internal::PbiIndexedBamReaderPrivate> d_;
};

} // namespace internal
} // namespace BAM

#endif // PBIINDEXEDBAMREADER_H
