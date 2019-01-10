// File Description
/// \file PbiIndexedBamReader.cpp
/// \brief Implements the PbiIndexedBamReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/PbiIndexedBamReader.h"

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <htslib/bgzf.h>
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include "MemoryUtils.h"
#include "PbiIndexIO.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

class PbiIndexedBamReader::PbiIndexedBamReaderPrivate
{
public:
    explicit PbiIndexedBamReaderPrivate(const std::string& pbiFilename)
        : pbiFilename_{pbiFilename}, currentBlockReadCount_{0}, numMatchingReads_{0}
    {
        PbiIndexIO io{pbiFilename, {PbiFile::Field::VIRTUAL_OFFSET}};
        header_ = io.Header();
    }

    void Filter(const PbiFilter newFilter)
    {
        // store request & reset counters
        filter_ = std::move(newFilter);
        currentBlockReadCount_ = 0;
        blocks_.clear();
        numMatchingReads_ = 0;

        // empty PBI (no reads)
        if (header_.numReads == 0) return;

        // maybe load/reload index
        const auto currentFilterFields = filter_.RequiredFields();
        const auto newFilterFields = newFilter.RequiredFields();
        const bool shouldLoadFields = !std::equal(
            currentFilterFields.cbegin(), currentFilterFields.cend(), newFilterFields.cbegin());
        if (!index_.is_initialized() || shouldLoadFields) {
            PbiIndexIO io{pbiFilename_, newFilterFields};
            index_ = io.Load();
        }

        // empty filter (all reads)
        if (filter_.IsEmpty()) {
            numMatchingReads_ = header_.numReads;
            blocks_.emplace_back(0, header_.numReads);
        }

        // apply filter, store contiguous blocks of reads
        else {

            IndexList indices;
            indices.reserve(header_.numReads);
            const auto& index = index_.get();
            for (size_t i = 0; i < header_.numReads; ++i) {
                if (filter_.Accepts(index, i)) {
                    indices.push_back(i);
                    ++numMatchingReads_;
                }
            }
            blocks_ = MergedIndexBlocks(std::move(indices));
        }

        // add virtual offsets to blocks
        const auto& fileOffsets = index_->BasicData().fileOffset_;
        for (IndexResultBlock& block : blocks_)
            block.virtualOffset_ = fileOffsets.at(block.firstIndex_);
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

    std::string pbiFilename_;
    PbiFilter filter_;
    PbiHeader header_;
    boost::optional<PbiRawData> index_;  //
    IndexResultBlocks blocks_;
    size_t currentBlockReadCount_;
    uint32_t numMatchingReads_;
};

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
    , d_{std::make_unique<PbiIndexedBamReaderPrivate>(File().PacBioIndexFilename())}
{
}

PbiIndexedBamReader::~PbiIndexedBamReader() = default;

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
