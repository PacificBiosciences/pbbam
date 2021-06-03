#ifndef PBBAM_KSEQREADER_H
#define PBBAM_KSEQREADER_H

#include <pbbam/Config.h>

#include <memory>
#include <string>

#include <htslib/kseq.h>
#include <zlib.h>

#include <pbbam/Deleters.h>

namespace PacBio {
namespace BAM {

class KSeqReader
{
public:
    explicit KSeqReader(const std::string& fn);

    KSeqReader(const KSeqReader&) = delete;
    KSeqReader(KSeqReader&&) noexcept;
    KSeqReader& operator=(const KSeqReader&) = delete;
    KSeqReader& operator=(KSeqReader&&) noexcept;
    virtual ~KSeqReader();

    std::string Name() const;
    std::string Bases() const;
    std::string Qualities() const;

    bool ReadNext();

private:
    // specialize kseq for gzstream (handles all of our types)
    KSEQ_INIT(gzFile, gzread)
    struct KSeqDeleter
    {
        void operator()(kseq_t* seq) const noexcept
        {
            if (seq) {
                kseq_destroy(seq);
            }
            seq = nullptr;
        }
    };

    std::unique_ptr<gzFile_s, GzFileDeleter> fp_;
    std::unique_ptr<kseq_t, KSeqDeleter> seq_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_KSEQREADER_H
