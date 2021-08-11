#ifndef PBBAM_ZMWCHUNKEDFASTXBGZFREADER_H
#define PBBAM_ZMWCHUNKEDFASTXBGZFREADER_H

#include <pbbam/Config.h>

#include "ZmwChunkedFastxReaderImpl.h"

#include <memory>

#include <htslib/kseq.h>

#include <pbbam/Deleters.h>

namespace PacBio {
namespace BAM {

class ZmwChunkedFastxBgzfReader final : public ZmwChunkedFastxReaderImpl
{
public:
    ZmwChunkedFastxBgzfReader(std::string filename, const size_t numChunks);

    void Seek(uint64_t pos) final;
    FastaSequence ReadNextFasta(bool skipName) final;
    FastqSequence ReadNextFastq(bool skipName) final;

private:
    int FetchRecord(bool getName);

    // specialize kseq_t for BGZF handle
    KSEQ_INIT(BGZF*, bgzf_read);
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

    std::unique_ptr<BGZF, HtslibBgzfDeleter> file_;
    std::unique_ptr<kseq_t, KSeqDeleter> seq_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWCHUNKEDFASTXBGZFREADER_H
