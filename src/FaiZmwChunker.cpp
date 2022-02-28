#include "PbbamInternalConfig.h"

#include "FaiZmwChunker.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

#include <cassert>

namespace PacBio {
namespace BAM {
namespace {

int32_t HoleNumber(const std::string& name)
{
    const auto firstSlash = name.find('/');
    if (firstSlash == std::string::npos) {
        throw std::runtime_error{
            "[pbbam] FAI chunking ERROR: could not parse hole number from name: " + name};
    }

    auto numberEnd = name.find('/', firstSlash + 1);
    if (numberEnd == std::string::npos) {
        numberEnd = name.size();
    }

    return std::stoi(name.substr(firstSlash + 1, (numberEnd - firstSlash)));
}

}  // namespace

FaiZmwChunker::FaiZmwChunker(const FaiIndex& index, const size_t numChunks)
{
    // zero chunks is error
    if (numChunks == 0) {
        throw std::runtime_error{
            "[pbbam] FAI chunking ERROR: requested chunk count must be greater than zero"};
    }

    // empty index is not (?), but quick return
    const auto& names = index.Names();
    if (names.empty()) {
        return;
    }

    // tease apart unique ZMWs
    int32_t currentHoleNumber = -1;
    std::vector<FaiZmwChunk> rawChunks;
    for (const auto& name : names) {
        const int32_t holeNumber = HoleNumber(name);
        if (holeNumber != currentHoleNumber) {
            rawChunks.emplace_back(FaiZmwChunk{name, index.Entry(name).SeqOffset, 1, 1});
            currentHoleNumber = holeNumber;
        } else {
            ++rawChunks.back().NumRecords;
        }
    }

    // no empty chunks (e.g. reduce the requested number, if small ZMW input)
    const size_t actualNumChunks = std::min(numChunks, rawChunks.size());

    // determine how many ZMWs should land in each chunk, spread roughly evenly
    const int minimum = (rawChunks.size() / actualNumChunks);
    const int modulo = (rawChunks.size() % actualNumChunks);
    std::vector<size_t> chunkCounts(actualNumChunks, minimum);
    for (int i = 0; i < modulo; ++i) {
        ++chunkCounts.at(i);
    }

    // collate zmw data into larger chunks
    size_t begin = 0;
    size_t end = 0;
    for (const auto n : chunkCounts) {

        // shift end down for this chunk
        end += n;
        assert(end <= rawChunks.size());

        // add data for this chunk
        FaiZmwChunk result = rawChunks.at(begin);
        result.NumZmws = n;
        for (size_t j = begin + 1; j < end; ++j) {
            result.NumRecords += rawChunks.at(j).NumRecords;
        }
        chunks_.emplace_back(std::move(result));

        // slide to next chunk
        begin = end;
    }
}

FaiZmwChunker::FaiZmwChunker(const std::string& filename, const size_t numChunks)
    : FaiZmwChunker{FaiIndex{filename}, numChunks}
{}

FaiZmwChunker::FaiZmwChunker(const FaiZmwChunker&) = default;

FaiZmwChunker::FaiZmwChunker(FaiZmwChunker&&) noexcept = default;

FaiZmwChunker& FaiZmwChunker::operator=(const FaiZmwChunker&) = default;

FaiZmwChunker& FaiZmwChunker::operator=(FaiZmwChunker&&) noexcept = default;

FaiZmwChunker::~FaiZmwChunker() = default;

const FaiZmwChunk& FaiZmwChunker::Chunk(size_t chunk) const { return chunks_.at(chunk); }

size_t FaiZmwChunker::NumChunks() const { return chunks_.size(); }

}  // namespace BAM
}  // namespace PacBio
