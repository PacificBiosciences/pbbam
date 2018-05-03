// File Description
/// \file FastqReader.cpp
/// \brief Implements the FastqReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqReader.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <htslib/faidx.h>

#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {
namespace internal {

struct FastqReaderPrivate
{
public:
    explicit FastqReaderPrivate(const std::string& fn) : stream_{fn}
    {
        if (!stream_)
            throw std::runtime_error{"FastqReader - could not open " + fn + " for reading"};
        FetchNext();
    }

    bool GetNext(FastqSequence& record)
    {
        if (name_.empty() && bases_.empty() && quals_.empty()) return false;
        record = FastqSequence{name_, bases_, quals_};
        FetchNext();
        return true;
    }

private:
    void FetchNext()
    {
        name_.clear();
        bases_.clear();
        quals_.clear();

        if (!stream_ || stream_.eof()) return;

        SkipNewlines();

        ReadName();
        ReadBases();
        stream_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // ignore "comment line"
        ReadQuals();
    }

    inline void SkipNewlines()
    {
        if (stream_.peek() == '\n')
            stream_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    void ReadName()
    {
        if (stream_.get() == '@') std::getline(stream_, name_, '\n');
    }

    void ReadBases() { std::getline(stream_, bases_, '\n'); }

    void ReadQuals() { std::getline(stream_, quals_, '\n'); }

private:
    std::ifstream stream_;
    std::string name_;
    std::string bases_;
    std::string quals_;
};

}  // namespace internal

FastqReader::FastqReader(const std::string& fn)
    : d_{std::make_unique<internal::FastqReaderPrivate>(fn)}
{
}

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
