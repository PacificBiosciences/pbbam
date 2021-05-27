#ifndef PBBAM_DELETERS_H
#define PBBAM_DELETERS_H

#include <pbbam/Config.h>

#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <zlib.h>

namespace PacBio {
namespace BAM {

// intended for use with std::shared_ptr<T>, std::unique_ptr<T>, etc

struct GzFileDeleter
{
    void operator()(gzFile fp) const noexcept
    {
        if (fp) {
            gzclose(fp);
        }
        fp = nullptr;
    }
};

struct HtslibBgzfDeleter
{
    void operator()(BGZF* bgzf) const noexcept
    {
        bgzf_close(bgzf);
        bgzf = nullptr;
    }
};

struct HtslibFastaIndexDeleter
{
    void operator()(faidx_t* fai) const noexcept
    {
        fai_destroy(fai);
        fai = nullptr;
    }
};

struct HtslibFileDeleter
{
    void operator()(samFile* file) const noexcept
    {
        if (file) {
            sam_close(file);
        }
        file = nullptr;
    }
};

struct HtslibHeaderDeleter
{
    void operator()(bam_hdr_t* hdr) const noexcept
    {
        bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

struct HtslibIndexDeleter
{
    void operator()(hts_idx_t* index) const noexcept
    {
        hts_idx_destroy(index);
        index = nullptr;
    }
};

struct HtslibIteratorDeleter
{
    void operator()(hts_itr_t* iter) const noexcept
    {
        hts_itr_destroy(iter);
        iter = nullptr;
    }
};

struct HtslibRecordDeleter
{
    void operator()(bam1_t* b) const noexcept
    {
        bam_destroy1(b);
        b = nullptr;
    }
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_DELETERS_H
