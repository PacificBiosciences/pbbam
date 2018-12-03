// File Description
/// \file FastqReader.cpp
/// \brief Implements the FastqReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqReader.h"

#include <memory>
#include <stdexcept>
#include <string>

#include <htslib/kseq.h>
#include <zlib.h>

#include "pbbam/FastqSequence.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

class FastqReader::FastqReaderPrivate
{
    KSEQ_INIT(gzFile, gzread)
    struct KSeqDeleter
    {
        void operator()(kseq_t* seq)
        {
            if (seq) kseq_destroy(seq);
            seq = nullptr;
        }
    };

public:
    explicit FastqReaderPrivate(const std::string& fn)
        : fp_{gzopen(fn.c_str(), "r")}, seq_{kseq_init(fp_)}
    {
        if (fp_ == nullptr || seq_.get() == nullptr)
            throw std::runtime_error{"Could not open " + fn + " for reading"};
    }

    ~FastqReaderPrivate() { gzclose(fp_); }

    bool GetNext(FastqSequence& record)
    {
        const auto result = kseq_read(seq_.get());
        if (result == -1)  // EOF
            return false;
        record = FastqSequence{std::string{seq_->name.s, seq_->name.l},
                               std::string{seq_->seq.s, seq_->seq.l},
                               std::string{seq_->qual.s, seq_->qual.l}};
        return true;
    }

private:
    gzFile fp_;
    std::unique_ptr<kseq_t, KSeqDeleter> seq_;
};

FastqReader::FastqReader(const std::string& fn) : d_{std::make_unique<FastqReaderPrivate>(fn)} {}

FastqReader::~FastqReader() {}

bool FastqReader::GetNext(FastqSequence& record) { return d_->GetNext(record); }

std::vector<FastqSequence> FastqReader::ReadAll(const std::string& fn)
{
    std::vector<FastqSequence> result;
    result.reserve(256);
    FastqReader reader{fn};
    FastqSequence s;
    while (reader.GetNext(s))
        result.emplace_back(s);
    return result;
}

}  // namespace BAM
}  // namespace PacBio
