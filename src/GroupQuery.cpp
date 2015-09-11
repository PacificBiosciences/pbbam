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

// Author: Yuan Li

#include "pbbam/GroupQuery.h"
#include "MemoryUtils.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

SequentialGroupQueryBase::SequentialGroupQueryBase(const BamFile & file) 
    : GroupQueryBase(file)
    , htsFile_(nullptr)
    , htsHeader_(nullptr)
    , nextRecord_()
{
    htsFile_.reset(sam_open(file.Filename().c_str(), "rb"), internal::HtslibFileDeleter());
    if (!htsFile_) 
        throw std::runtime_error("could not open BAM file for reading");

    htsHeader_.reset(sam_hdr_read(htsFile_.get()), internal::HtslibHeaderDeleter());
    if (!htsHeader_) 
        throw std::runtime_error("could not read BAM header data");
}

bool SequentialGroupQueryBase::GetNext(vector<BamRecord> & records) 
{
    records.clear();

    if (nextRecord_.Impl().Name() != "") {
        records.push_back(nextRecord_);
        nextRecord_ = BamRecord();
    }

    while(true) {
        BamRecord record(file_.Header());
        const int result = sam_read1(htsFile_.get(),
                                     htsHeader_.get(),
                                     internal::BamRecordMemory::GetRawData(record).get());
        internal::BamRecordMemory::UpdateRecordTags(record);
        if (result >= 0) { // get next record
            if (records.size() == 0) {
                records.push_back(record); // add the first record
            } else {
                if (InSameGroup(record, records[0])) {
                    records.push_back(record); // add remaining record
                } else {
                    nextRecord_ = record; // store record from another zmw
                    return true;
                }
            }
        } else { // unable to get next record
            if (records.size() > 0) return true; // Has records to return
            else return false; // Has no records to return
        }
    }
    assert(false); // Should not reach here.
    return false;
}
