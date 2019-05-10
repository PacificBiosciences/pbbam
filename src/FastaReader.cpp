// File Description
/// \file FastaReader.cpp
/// \brief Implements the FastaReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastaReader.h"

#include <memory>
#include <stdexcept>
#include <string>

#include <htslib/kseq.h>
#include <zlib.h>

#include "pbbam/FastaSequence.h"

namespace PacBio {
namespace BAM {

class FastaReader::FastaReaderPrivate
{
    KSEQ_INIT(gzFile, gzread)
    struct KSeqDeleter
    {
        void operator()(kseq_t* seq) const
        {
            if (seq) kseq_destroy(seq);
            seq = nullptr;
        }
    };

public:
    explicit FastaReaderPrivate(const std::string& fn)
        : fp_{gzopen(fn.c_str(), "r")}, seq_{kseq_init(fp_)}
    {
        if (fp_ == nullptr || seq_.get() == nullptr)
            throw std::runtime_error{"FastaReader: could not open file for reading: " + fn};
    }

    ~FastaReaderPrivate() { gzclose(fp_); }

    bool GetNext(FastaSequence& record)
    {
        const auto result = kseq_read(seq_.get());
        if (result == -1)  // EOF
            return false;
        record = FastaSequence{std::string{seq_->name.s, seq_->name.l},
                               std::string{seq_->seq.s, seq_->seq.l}};
        return true;
    }

private:
    gzFile fp_;
    std::unique_ptr<kseq_t, KSeqDeleter> seq_;
};

FastaReader::FastaReader(const std::string& fn)
    : internal::QueryBase<FastaSequence>{}, d_{std::make_unique<FastaReaderPrivate>(fn)}
{
}

FastaReader::FastaReader(FastaReader&&) = default;

FastaReader& FastaReader::operator=(FastaReader&&) = default;

FastaReader::~FastaReader() = default;

bool FastaReader::GetNext(FastaSequence& record) { return d_->GetNext(record); }

std::vector<FastaSequence> FastaReader::ReadAll(const std::string& fn)
{
    std::vector<FastaSequence> result;
    result.reserve(256);
    FastaReader reader{fn};
    FastaSequence s;
    while (reader.GetNext(s))
        result.emplace_back(s);
    return result;
}

}  // namespace BAM
}  // namespace PacBio
