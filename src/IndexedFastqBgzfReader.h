#ifndef PBBAM_INDEXEDFASTQBGZFREADER_H
#define PBBAM_INDEXEDFASTQBGZFREADER_H

#include <pbbam/Config.h>

#include "IndexedFastqReaderImpl.h"

#include <memory>

#include <htslib/kseq.h>
#include <pbcopper/data/Position.h>
#include <pbcopper/data/QualityValues.h>

#include <pbbam/Deleters.h>

namespace PacBio {
namespace BAM {

class IndexedFastqBgzfReader final : public IndexedFastqReaderImpl
{
public:
    IndexedFastqBgzfReader(std::string filename);

    std::pair<std::string, Data::QualityValues> Subsequence(const std::string& id,
                                                            Data::Position start,
                                                            Data::Position end) final;

private:
    int FetchRecord();

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

#endif  // PBBAM_INDEXEDFASTQBGZFREADER_H
