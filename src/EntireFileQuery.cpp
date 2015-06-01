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

#include "pbbam/EntireFileQuery.h"
#include "pbbam/BamFile.h"

#include "pbbam/internal/SequentialMergeStrategy.h"

#include "MemoryUtils.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

EntireFileQuery::EntireFileQuery(const BamFile& file)
    : QueryBase(file)
    , htsFile_(nullptr)
    , htsHeader_(nullptr)
{
    htsFile_.reset(sam_open(file.Filename().c_str(), "rb"), internal::HtslibFileDeleter());
    if (!htsFile_)
        throw std::exception();

    htsHeader_.reset(sam_hdr_read(htsFile_.get()), internal::HtslibHeaderDeleter());
    if (!htsHeader_)
        throw std::exception();
}

bool EntireFileQuery::GetNext(BamRecord& record)
{
    const int result = sam_read1(htsFile_.get(),
                                 htsHeader_.get(),
                                 internal::BamRecordMemory::GetRawData(record).get());
    // success
    if (result >= 0)
        return true;

    // normal EOF
    else if (result == -1)
        return false;

    // error (truncated file, etc)
    else
        throw std::exception();
}

namespace PacBio {
namespace BAM {
namespace staging {

class EntireFileIterator : public internal::IBamFileIterator
{
public:
    EntireFileIterator(const BamFile& bamFile)
        : internal::IBamFileIterator(bamFile)
    {
        htsFile_.reset(sam_open(bamFile.Filename().c_str(), "rb"));
        if (!htsFile_)
            throw std::exception();

        htsHeader_.reset(sam_hdr_read(htsFile_.get()));
        if (!htsHeader_)
            throw std::exception();
    }

public:
    bool GetNext(BamRecord& record) {

        const int result = sam_read1(htsFile_.get(),
                                     htsHeader_.get(),
                                     internal::BamRecordMemory::GetRawData(record).get());
        // success
        if (result >= 0)
            return true;

        // normal EOF
        else if (result == -1)
            return false;

        // error (truncated file, etc)
        else
            throw std::exception();
    }

private:
    unique_ptr<samFile,   internal::HtslibFileDeleter>   htsFile_;
    unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> htsHeader_;
};

EntireFileQuery::EntireFileQuery(const DataSet& dataset)
    : internal::IQuery(dataset)
{
    // check files
    // if SO all coordinate
    // else if SO all queryname
    // else SO unsorted/unknown
    mergeStrategy_.reset(new internal::SequentialMergeStrategy(CreateIterators()));
}

EntireFileQuery::FileIterPtr EntireFileQuery::CreateIterator(const BamFile& bamFile)
{ return FileIterPtr(new EntireFileIterator(bamFile)); }

} // namespace staging
} // namespace BAM
} // namspace PacBio
