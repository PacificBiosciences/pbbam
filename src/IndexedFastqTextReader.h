#ifndef PBBAM_INDEXEDFASTQTEXTREADER_H
#define PBBAM_INDEXEDFASTQTEXTREADER_H

#include <pbbam/Config.h>

#include "IndexedFastqReaderImpl.h"

#include <pbcopper/utility/Deleters.h>

#include <htslib/kseq.h>

#include <memory>

#include <cstdio>

namespace PacBio {
namespace BAM {

class IndexedFastqTextReader final : public IndexedFastqReaderImpl
{
public:
    IndexedFastqTextReader(std::string filename);

    std::pair<std::string, Data::QualityValues> Subsequence(const std::string& id,
                                                            Data::Position start,
                                                            Data::Position end) final;

private:
    int FetchRecord();

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

#endif  // PBBAM_INDEXEDFASTQTEXTREADER_H
