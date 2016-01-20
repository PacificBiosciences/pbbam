// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#ifndef UNMAPPEDREADSQUERY_H
#define UNMAPPEDREADSQUERY_H

// TODO: Not yet implemented.
//       Can't get bam_itr_queryi(idx, HTS_IDX_NOCOOR, -1, -1) to work reliably at the moment.

//#include "pbbam/QueryBase.h"
//#include <string>

//namespace PacBio {
//namespace BAM {

//class BamFile;

//class PBBAM_EXPORT UnmappedReadsQuery : public QueryBase
//{
//public:
//    UnmappedReadsQuery(const BamFile& file);

//protected:
//    bool GetNext(BamRecord& record);

//private:
//    std::shared_ptr<samFile>   file_;
//    std::shared_ptr<bam_hdr_t> header_;
//    std::shared_ptr<hts_idx_t> index_;
//    std::shared_ptr<hts_itr_t> iterator_;
//};

//} // namespace BAM
//} // namspace PacBio

#endif // UNMAPPEDREADSQUERY_H
