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

#ifndef MEMORYUTILS_H
#define MEMORYUTILS_H

#include "pbbam/Config.h"
#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/BamRecordImpl.h"
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <memory>

namespace PacBio {
namespace BAM {

class BamHeader;

namespace internal {

// intended for use with PBBAM_SHARED_PTR<T>, std::unique_ptr<T>, etc

struct HtslibBgzfDeleter
{
    void operator()(BGZF* bgzf)
    {
        if (bgzf) 
            bgzf_close(bgzf);
        bgzf = nullptr;
    }
};

struct HtslibFileDeleter
{
    void operator()(samFile* file)
    {
        if (file)
            sam_close(file);
        file = nullptr;
    }
};

struct HtslibHeaderDeleter
{
    void operator()(bam_hdr_t* hdr)
    {
        if (hdr)
            bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

struct HtslibIndexDeleter
{
    void operator()(hts_idx_t* index)
    {
        if (index)
            hts_idx_destroy(index);
        index = nullptr;
    }
};

struct HtslibIteratorDeleter
{
    void operator()(hts_itr_t* iter)
    {
        if (iter)
            hts_itr_destroy(iter);
        iter = nullptr;
    }
};

struct HtslibRecordDeleter
{
    void operator()(bam1_t* b)
    {
        if (b)
            bam_destroy1(b);
        b = nullptr;
    }
};

class BamHeaderMemory
{
public:
    static BamHeader FromRawData(bam_hdr_t* header);
    static PBBAM_SHARED_PTR<bam_hdr_t> MakeRawHeader(const BamHeader& header);
//    static PBBAM_SHARED_PTR<bam_hdr_t> MakeRawHeader(const BamHeader& header);
};

class BamRecordMemory
{
public:
    static const BamRecordImpl& GetImpl(const BamRecord& r);
    static const BamRecordImpl& GetImpl(const BamRecord* r);
    static PBBAM_SHARED_PTR<bam1_t> GetRawData(const BamRecord& r);
    static PBBAM_SHARED_PTR<bam1_t> GetRawData(const BamRecord* r);
    static PBBAM_SHARED_PTR<bam1_t> GetRawData(const BamRecordImpl& impl);
    static PBBAM_SHARED_PTR<bam1_t> GetRawData(const BamRecordImpl* impl);

    static void UpdateRecordTags(const BamRecord& r);
    static void UpdateRecordTags(const BamRecordImpl& r);
};

inline const BamRecordImpl& BamRecordMemory::GetImpl(const BamRecord& r)
{ return r.impl_; }

inline const BamRecordImpl& BamRecordMemory::GetImpl(const BamRecord* r)
{ return r->impl_; }

inline PBBAM_SHARED_PTR<bam1_t> BamRecordMemory::GetRawData(const BamRecord& r)
{ return GetRawData(r.impl_); }

inline PBBAM_SHARED_PTR<bam1_t> BamRecordMemory::GetRawData(const BamRecord* r)
{ return GetRawData(r->impl_); }

inline PBBAM_SHARED_PTR<bam1_t> BamRecordMemory::GetRawData(const BamRecordImpl& impl)
{ return impl.d_; }

inline PBBAM_SHARED_PTR<bam1_t> BamRecordMemory::GetRawData(const BamRecordImpl* impl)
{ return impl->d_; }

inline void BamRecordMemory::UpdateRecordTags(const BamRecord& r)
{ UpdateRecordTags(r.impl_); }

inline void BamRecordMemory::UpdateRecordTags(const BamRecordImpl& r)
{ r.UpdateTagMap(); }

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // MEMORYUTILS_H
