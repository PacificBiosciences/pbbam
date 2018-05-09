// File Description
/// \file FastaReader.cpp
/// \brief Implements the FastaReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastaReader.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <htslib/faidx.h>

#include "pbbam/MakeUnique.h"
#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {
namespace internal {

struct FastaReaderPrivate
{
    std::ifstream stream_;
    std::string name_;
    std::string bases_;

    FastaReaderPrivate(const std::string& fn) : stream_{fn}
    {
        if (!stream_)
            throw std::runtime_error{"FastaReader - could not open " + fn + " for reading"};
        FetchNext();
    }

    bool GetNext(FastaSequence& record)
    {
        if (name_.empty() && bases_.empty()) return false;
        record = FastaSequence{name_, bases_};
        FetchNext();
        return true;
    }

private:
    void FetchNext()
    {
        name_.clear();
        bases_.clear();

        SkipNewlines();
        ReadName();
        ReadBases();

        bases_ = RemoveAllWhitespace(std::move(bases_));
    }

    inline void SkipNewlines()
    {
        if (!stream_) return;
        if (stream_.peek() == '\n')
            stream_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    void ReadName()
    {
        if (!stream_) return;
        if (stream_.get() == '>') std::getline(stream_, name_, '\n');
    }

    void ReadBases()
    {
        if (!stream_) return;
        int p = stream_.peek();
        while (static_cast<char>(p) != '>' && p != EOF) {
            if (!stream_) return;
            std::string line;
            std::getline(stream_, line, '\n');
            bases_ += line;
            if (!stream_) return;
            p = stream_.peek();
        }
    }
};

}  // namespace internal

FastaReader::FastaReader(const std::string& fn)
    : d_{std::make_unique<internal::FastaReaderPrivate>(fn)}
{
}

FastaReader::~FastaReader() {}

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
