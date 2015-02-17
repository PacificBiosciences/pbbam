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

//#include "pbbam/UnmappedReadsQuery.h"
//#include "pbbam/BamFile.h"
//#include "MemoryUtils.h"

//#include <iostream>

//using namespace PacBio;
//using namespace PacBio::BAM;
//using namespace std;

//UnmappedReadsQuery::UnmappedReadsQuery(const BamFile& file)
//    : QueryBase()
//{
//    // open file
//    file_.reset(sam_open(file.Filename().c_str(), "rb"), internal::HtslibFileDeleter());
//    if (!file_) {
//        error_ = UnmappedReadsQuery::FileOpenError;
//        return;
//    }

//    // open index
//    index_.reset(bam_index_load(file.Filename().c_str()), internal::HtslibIndexDeleter());
//    if (!index_) {
//        error_ = UnmappedReadsQuery::IndexFileOpenError;
//        return;
//    }

//    // initialize query
//    iterator_.reset(bam_itr_queryi(index_.get(), HTS_IDX_NOCOOR, 0, 0), internal::HtslibIteratorDeleter());
//    if (iterator_) {

//        cerr << endl
//             << "UnmappedQueryReads::iterator" << endl
//             << "read_rest: " << iterator_->read_rest << endl
//             << "finished: " << iterator_->finished << endl
//             << "dummy: " << iterator_->dummy << endl
//             << "tid: " << iterator_->tid << endl
//             << "beg: " << iterator_->beg << endl
//             << "end: " << iterator_->end << endl
//             << "n_off: " << iterator_->n_off << endl
//             << "i: " << iterator_->i << endl
//             << "curr_off: " << iterator_->curr_off << endl
//             << endl;


////        uint32_t read_rest:1, finished:1, dummy:29;
////        int tid, beg, end, n_off, i;
////        uint64_t curr_off;
////        hts_pair64_t *off;
////        hts_readrec_func *readrec;
////        struct {
////            int n, m;
////            int *a;
////        } bins;

//    }
//}

//bool UnmappedReadsQuery::GetNext(BamRecord& record)
//{
//    if (error_ == UnmappedReadsQuery::NoError && iterator_) {
//        const int result = bam_itr_next(file_.get(), iterator_.get(), record.RawData().get());
//        if ( result > 0 )
//            return true;
//        else {
//            cerr << "ERROR - result: " << result << endl;
//            if ( result == -4 ) {

//                bam1_t* b = record.RawData().get();
//                bam1_core_t* c = &b->core;
//                bool nonBgzfErrorFound = false;

//                if (b->l_data < 0) {
//                    cerr << "ERROR: bam1_t::l_data < 0" << endl;
//                    nonBgzfErrorFound = true;
//                }
//                if (c->l_qseq < 0) {
//                    cerr << "ERROR: bam1_t::core::l_qseq < 0" << endl;
//                    nonBgzfErrorFound = true;
//                }
//                if (!b->data) {
//                    cerr << "ERROR: bam1_t::data is null" << endl;
//                    nonBgzfErrorFound = true;
//                }
//                if  (!nonBgzfErrorFound)
//                    cerr << "ERROR: in bam_read1(), bgzf_read(fp, b->data, b->l_data) returned unexpected value" << endl;
//            }
//        }
//    }
//    else {
//        cerr << "UnmappedReadsQuery::HasError() - " << Error() << endl;
//    }


//    return false;
//}
