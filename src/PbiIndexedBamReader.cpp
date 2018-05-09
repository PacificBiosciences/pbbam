// File Description
/// \file PbiIndexedBamReader.cpp
/// \brief Implements the PbiIndexedBamReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/PbiIndexedBamReader.h"

#include <cstddef>
#include <cstdint>
#include <iostream>

#include <htslib/bgzf.h>

#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {
namespace internal {

struct PbiIndexedBamReaderPrivate
{
public:
    PbiIndexedBamReaderPrivate(const std::string& pbiFilename)
        : index_{pbiFilename}, currentBlockReadCount_{0}, numMatchingReads_{0}
    {
    }

    void ApplyOffsets()
    {
        const auto& fileOffsets = index_.BasicData().fileOffset_;
        for (IndexResultBlock& block : blocks_)
            block.virtualOffset_ = fileOffsets.at(block.firstIndex_);
    }

    void Filter(const PbiFilter filter)
    {
        // store request & reset counters
        filter_ = std::move(filter);
        currentBlockReadCount_ = 0;
        blocks_.clear();
        numMatchingReads_ = 0;

        // find blocks of reads passing filter criteria
        const auto totalReads = index_.NumReads();
        if (totalReads == 0) {  // empty PBI - no reads to use
            return;
        } else if (filter_.IsEmpty()) {  // empty filter - use all reads
            numMatchingReads_ = totalReads;
            blocks_.emplace_back(0, totalReads);
        } else {
            IndexList indices;
            indices.reserve(totalReads);
            for (size_t i = 0; i < totalReads; ++i) {
                if (filter_.Accepts(index_, i)) {
                    indices.push_back(i);
                    ++numMatchingReads_;
                }
            }
            blocks_ = MergedIndexBlocks(std::move(indices));
        }

        // apply offsets
        ApplyOffsets();
    }

    IndexResultBlocks MergedIndexBlocks(IndexList indices) const
    {
        if (indices.empty()) return {};

        std::sort(indices.begin(), indices.end());
        auto newEndIter = std::unique(indices.begin(), indices.end());
        auto numIndices = std::distance(indices.begin(), newEndIter);
        auto result = IndexResultBlocks{IndexResultBlock{indices.at(0), 1}};
        for (auto i = 1; i < numIndices; ++i) {
            if (indices.at(i) == indices.at(i - 1) + 1)
                ++result.back().numReads_;
            else
                result.emplace_back(indices.at(i), 1);
        }
        return result;
    }

    int ReadRawData(BGZF* bgzf, bam1_t* b)
    {
        // no data to fetch, return false
        if (blocks_.empty()) return -1;  // "EOF"

        // if on new block, seek to its first record
        if (currentBlockReadCount_ == 0) {
            const auto seekResult = bgzf_seek(bgzf, blocks_.at(0).virtualOffset_, SEEK_SET);
            if (seekResult == -1) throw std::runtime_error{"could not seek in BAM file"};
        }

        // read next record
        const auto result = bam_read1(bgzf, b);

        // update counters. if block finished, pop & reset
        ++currentBlockReadCount_;
        if (currentBlockReadCount_ == blocks_.at(0).numReads_) {
            blocks_.pop_front();
            currentBlockReadCount_ = 0;
        }

        return result;
    }

public:
    PbiFilter filter_;
    PbiRawData index_;
    IndexResultBlocks blocks_;
    size_t currentBlockReadCount_;
    uint32_t numMatchingReads_;
};

}  // namespace internal

PbiIndexedBamReader::PbiIndexedBamReader(PbiFilter filter, const std::string& filename)
    : PbiIndexedBamReader{std::move(filter), BamFile{filename}}
{
}

PbiIndexedBamReader::PbiIndexedBamReader(PbiFilter filter, BamFile bamFile)
    : PbiIndexedBamReader{std::move(bamFile)}
{
    Filter(std::move(filter));
}

PbiIndexedBamReader::PbiIndexedBamReader(const std::string& bamFilename)
    : PbiIndexedBamReader{BamFile{bamFilename}}
{
}

PbiIndexedBamReader::PbiIndexedBamReader(BamFile bamFile)
    : BamReader{std::move(bamFile)}
    , d_{std::make_unique<internal::PbiIndexedBamReaderPrivate>(File().PacBioIndexFilename())}
{
}

PbiIndexedBamReader::~PbiIndexedBamReader() {}

int PbiIndexedBamReader::ReadRawData(BGZF* bgzf, bam1_t* b)
{
    assert(d_);
    return d_->ReadRawData(bgzf, b);
}

const PbiFilter& PbiIndexedBamReader::Filter() const
{
    assert(d_);
    return d_->filter_;
}

PbiIndexedBamReader& PbiIndexedBamReader::Filter(PbiFilter filter)
{
    assert(d_);
    d_->Filter(std::move(filter));
    return *this;
}

uint32_t PbiIndexedBamReader::NumReads() const
{
    assert(d_);
    return d_->numMatchingReads_;
}

}  // namespace BAM
}  // namespace PacBio
