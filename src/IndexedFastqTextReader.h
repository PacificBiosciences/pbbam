// File Description
/// \file IndexedFastqTextReader.h
/// \brief Defines the IndexedFastqTextReader class.
//
// Author: Derek Barnett

#ifndef INDEXEDFASTQTEXTREADER_H
#define INDEXEDFASTQTEXTREADER_H

#include "pbbam/Config.h"

#include "IndexedFastqReaderImpl.h"

#include <cstdio>

#include <memory>

#include <htslib/kseq.h>
#include <pbcopper/utility/Deleters.h>

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
            if (seq) kseq_destroy(seq);
            seq = nullptr;
        }
    };

    std::unique_ptr<FILE, Utility::FileDeleter> file_;
    std::unique_ptr<kseq_t, KSeqDeleter> seq_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // INDEXEDFASTQTEXTREADER_H
