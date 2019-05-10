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

#include "KSeqReader.h"

namespace PacBio {
namespace BAM {

class FastaReader::FastaReaderPrivate : public KSeqReader
{
public:
    explicit FastaReaderPrivate(const std::string& fn) : KSeqReader{fn} {}

    bool GetNext(FastaSequence& record)
    {
        const auto readOk = ReadNext();
        if (!readOk) return false;

        record = FastaSequence{std::string{seq_->name.s, seq_->name.l},
                               std::string{seq_->seq.s, seq_->seq.l}};
        return true;
    }
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
