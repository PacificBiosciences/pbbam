// File Description
/// \file FastqReader.cpp
/// \brief Implements the FastqReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqReader.h"

#include <stdexcept>

#include "KSeqReader.h"

namespace PacBio {
namespace BAM {

class FastqReader::FastqReaderPrivate : public KSeqReader
{
public:
    explicit FastqReaderPrivate(const std::string& fn) : KSeqReader{fn} {}

    bool GetNext(FastqSequence& record)
    {
        const auto readOk = ReadNext();
        if (!readOk) return false;

        record = FastqSequence{std::string{seq_->name.s, seq_->name.l},
                               std::string{seq_->seq.s, seq_->seq.l},
                               std::string{seq_->qual.s, seq_->qual.l}};
        return true;
    }
};

FastqReader::FastqReader(const std::string& fn)
    : internal::QueryBase<FastqSequence>{}, d_{std::make_unique<FastqReaderPrivate>(fn)}
{
}

FastqReader::FastqReader(FastqReader&&) = default;

FastqReader& FastqReader::operator=(FastqReader&&) = default;

FastqReader::~FastqReader() = default;

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
