#ifndef PBBAM_ZMWCHUNKEDFASTXTEXTREADER_H
#define PBBAM_ZMWCHUNKEDFASTXTEXTREADER_H

#include <pbbam/Config.h>

#include "ZmwChunkedFastxReaderImpl.h"

#include <pbcopper/utility/Deleters.h>

#include <htslib/kseq.h>

#include <memory>

#include <cstdio>

namespace PacBio {
namespace BAM {

class ZmwChunkedFastxTextReader final : public ZmwChunkedFastxReaderImpl
{
public:
    ZmwChunkedFastxTextReader(std::string filename, size_t numChunks);

    void Seek(uint64_t pos) final;
    FastaSequence ReadNextFasta(bool skipName) final;
    FastqSequence ReadNextFastq(bool skipName) final;

private:
    int FetchRecord(bool getName);

    // kseq needs a '__read' function with this signature, so fread does not work
    // in this case. gzread/bgzf_read match but we want better seek performance
    // than gzstream and are specifically not using indexed BGZF
    static int ReadFromFile(FILE* fp, void* data, size_t length);

    // specialize kseq_t for FILE handle
    KSEQ_INIT(FILE*, ReadFromFile)
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

    std::unique_ptr<FILE, Utility::FileDeleter> file_;
    std::unique_ptr<kseq_t, KSeqDeleter> seq_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWCHUNKEDFASTXTEXTREADER_H
