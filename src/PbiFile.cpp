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

#include "pbbam/PbiFile.h"
#include "pbbam/BamFile.h"
#include "pbbam/PbiBuilder.h"
#include "MemoryUtils.h"
#include <htslib/sam.h>
#include <cassert>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::PbiFile;
using namespace std;

namespace PacBio {
namespace BAM {
namespace PbiFile {

void CreateFrom(const BamFile& bamFile)
{
    // open input file for file handle & header
    unique_ptr<samFile,internal::HtslibFileDeleter> htsFile(sam_open(bamFile.Filename().c_str(), "rb"));
    if (!htsFile)
        throw std::runtime_error("could not open BAM file for reading");

    unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> htsHeader(sam_hdr_read(htsFile.get()));
    if (!htsHeader)
        throw std::runtime_error("could not read BAM header data");

    samFile*   fp  = htsFile.get();
    bam_hdr_t* hdr = htsHeader.get();
    assert(fp);
    assert(hdr);

    // setup our record object
    BamRecord record;
    bam1_t* b = internal::BamRecordMemory::GetRawData(record).get();
    if (b == 0)
        throw std::runtime_error("could not allocate BAM record");

    // iterate through file, building index data
    PbiBuilder builder(bamFile.PacBioIndexFilename(), bamFile.Header().Sequences().size());
    int64_t offset = bgzf_tell(fp->fp.bgzf);
    int result = 0;
    while ((result = sam_read1(fp, hdr, b)) >= 0) {
        internal::BamRecordMemory::UpdateRecordTags(record);
        builder.AddRecord(record, offset);
        offset = bgzf_tell(fp->fp.bgzf);
    }
}

} // namespace PbiFile
} // namespace BAM
} // namespace PacBio
