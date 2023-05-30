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
    ZmwChunkedFastxTextReader(std::string filename, std::size_t numChunks);

    void Seek(std::uint64_t pos) final;
    FastaSequence ReadNextFasta(bool skipName) final;
    FastqSequence ReadNextFastq(bool skipName) final;

private:
    int FetchRecord(bool getName);

    // kseq needs a '__read' function with this signature, so std::fread does not work
    // in this case. gzread/bgzf_read match but we want better seek performance
    // than gzstream and are specifically not using indexed BGZF
    static int ReadFromFile(std::FILE* fp, void* data, std::size_t length);

    // specialize kseq_t for std::FILE handle
    KSEQ_INIT(std::FILE*, ReadFromFile)
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

    std::unique_ptr<std::FILE, Utility::FileDeleter> file_;
    std::unique_ptr<kseq_t, KSeqDeleter> seq_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWCHUNKEDFASTXTEXTREADER_H
