// Author: Derek Barnett

#ifndef MEMORYUTILS_H
#define MEMORYUTILS_H

#include <cstdio>
#include <memory>

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/BamRecordImpl.h"
#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

class BamHeader;

namespace internal {

// intended for use with std::shared_ptr<T>, std::unique_ptr<T>, etc

struct FileDeleter
{
    void operator()(std::FILE* fp)
    {
        if (fp) std::fclose(fp);
        fp = nullptr;
    }
};

struct HtslibBgzfDeleter
{
    void operator()(BGZF* bgzf)
    {
        if (bgzf) bgzf_close(bgzf);
        bgzf = nullptr;
    }
};

struct HtslibFileDeleter
{
    void operator()(samFile* file)
    {
        if (file) sam_close(file);
        file = nullptr;
    }
};

struct HtslibHeaderDeleter
{
    void operator()(bam_hdr_t* hdr)
    {
        if (hdr) bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

struct HtslibIndexDeleter
{
    void operator()(hts_idx_t* index)
    {
        if (index) hts_idx_destroy(index);
        index = nullptr;
    }
};

struct HtslibIteratorDeleter
{
    void operator()(hts_itr_t* iter)
    {
        if (iter) hts_itr_destroy(iter);
        iter = nullptr;
    }
};

struct HtslibRecordDeleter
{
    void operator()(bam1_t* b)
    {
        if (b) bam_destroy1(b);
        b = nullptr;
    }
};

class BamHeaderMemory
{
public:
    static BamHeader FromRawData(bam_hdr_t* header);
    static std::shared_ptr<bam_hdr_t> MakeRawHeader(const BamHeader& header);
};

class BamRecordMemory
{
public:
    static const BamRecordImpl& GetImpl(const BamRecord& r);
    static const BamRecordImpl& GetImpl(const BamRecord* r);
    static std::shared_ptr<bam1_t> GetRawData(const BamRecord& r);
    static std::shared_ptr<bam1_t> GetRawData(const BamRecord* r);
    static std::shared_ptr<bam1_t> GetRawData(const BamRecordImpl& impl);
    static std::shared_ptr<bam1_t> GetRawData(const BamRecordImpl* impl);

    static void UpdateRecordTags(const BamRecord& r);
    static void UpdateRecordTags(const BamRecordImpl& r);
};

inline const BamRecordImpl& BamRecordMemory::GetImpl(const BamRecord& r) { return r.impl_; }

inline const BamRecordImpl& BamRecordMemory::GetImpl(const BamRecord* r) { return r->impl_; }

inline std::shared_ptr<bam1_t> BamRecordMemory::GetRawData(const BamRecord& r)
{
    return GetRawData(r.impl_);
}

inline std::shared_ptr<bam1_t> BamRecordMemory::GetRawData(const BamRecord* r)
{
    return GetRawData(r->impl_);
}

inline std::shared_ptr<bam1_t> BamRecordMemory::GetRawData(const BamRecordImpl& impl)
{
    return impl.d_;
}

inline std::shared_ptr<bam1_t> BamRecordMemory::GetRawData(const BamRecordImpl* impl)
{
    return impl->d_;
}

inline void BamRecordMemory::UpdateRecordTags(const BamRecord& r) { UpdateRecordTags(r.impl_); }

inline void BamRecordMemory::UpdateRecordTags(const BamRecordImpl& r) { r.UpdateTagMap(); }

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // MEMORYUTILS_H
