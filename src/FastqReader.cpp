// File Description
/// \file FastqReader.cpp
/// \brief Implements the FastqReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqReader.h"

#include <stdexcept>

#include "FormatUtils.h"
#include "KSeqReader.h"

namespace PacBio {
namespace BAM {

class FastqReader::FastqReaderPrivate
{
public:
    explicit FastqReaderPrivate(const std::string& fn)
    {
        // validate extension
        if (!FormatUtils::IsFastqFilename(fn)) {
            throw std::runtime_error{"FastqReader: filename '" + fn +
                                     "' is not recognized as a FASTQ file."};
        }
        reader_ = std::make_unique<KSeqReader>(fn);
    }

    bool GetNext(FastqSequence& record)
    {
        const auto ok = reader_->ReadNext();
        if (!ok) return false;  // not error, could be EOF

        record = FastqSequence{reader_->Name(), reader_->Bases(), reader_->Qualities()};
        return true;
    }

    std::unique_ptr<KSeqReader> reader_;
};

FastqReader::FastqReader(const std::string& fn)
    : internal::QueryBase<FastqSequence>{}, d_{std::make_unique<FastqReaderPrivate>(fn)}
{
}

FastqReader::FastqReader(FastqReader&&) noexcept = default;

FastqReader& FastqReader::operator=(FastqReader&&) noexcept = default;

FastqReader::~FastqReader() = default;

bool FastqReader::GetNext(FastqSequence& record) { return d_->GetNext(record); }

std::vector<FastqSequence> FastqReader::ReadAll(const std::string& fn)
{
    std::vector<FastqSequence> result;
    result.reserve(256);
    FastqReader reader{fn};
    for (const auto& seq : reader)
        result.emplace_back(seq);
    return result;
}

}  // namespace BAM
}  // namespace PacBio
