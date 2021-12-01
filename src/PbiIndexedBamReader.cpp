#include "PbbamInternalConfig.h"

#include <pbbam/PbiIndexedBamReader.h>

#include <cstddef>
#include <cstdint>

#include <sstream>
#include <stdexcept>

#include <htslib/bgzf.h>

#include "ErrnoReason.h"

namespace PacBio {
namespace BAM {

class PbiIndexedBamReader::PbiIndexedBamReaderPrivate
{
public:
    explicit PbiIndexedBamReaderPrivate(BamFile file, const std::shared_ptr<PbiRawData>& index)
        : file_{std::move(file)}, index_{index}, currentBlockReadCount_{0}, numMatchingReads_{0}
    {
    }

    void ApplyOffsets()
    {
        const auto& fileOffsets = index_->BasicData().fileOffset_;
        for (IndexResultBlock& block : blocks_) {
            block.virtualOffset_ = fileOffsets.at(block.firstIndex_);
        }
    }

    void Filter(const PbiFilter filter)
    {
        // store request & reset counters
        filter_ = std::move(filter);
        currentBlockReadCount_ = 0;
        blocks_.clear();
        numMatchingReads_ = 0;

        // find blocks of reads passing filter criteria
        const auto totalReads = index_->NumReads();
        if (totalReads == 0) {  // empty PBI - no reads to use
            return;
        } else if (filter_.IsEmpty()) {  // empty filter - use all reads
            numMatchingReads_ = totalReads;
            blocks_.emplace_back(0, totalReads);
        } else {
            IndexList indices;
            indices.reserve(totalReads);
            const auto& idx = *index_;
            for (size_t i = 0; i < totalReads; ++i) {
                if (filter_.Accepts(idx, i)) {
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
        if (indices.empty()) {
            return {};
        }

        std::sort(indices.begin(), indices.end());
        const auto newEndIter = std::unique(indices.begin(), indices.end());
        const auto numIndices = std::distance(indices.begin(), newEndIter);
        auto result = IndexResultBlocks{IndexResultBlock{indices.at(0), 1}};
        for (auto i = 1; i < numIndices; ++i) {
            if (indices.at(i) == indices.at(i - 1) + 1) {
                ++result.back().numReads_;
            } else {
                result.emplace_back(indices.at(i), 1);
            }
        }
        return result;
    }

    int ReadRawData(BGZF* bgzf, bam1_t* b)
    {
        // no data to fetch, return false
        if (blocks_.empty()) {
            return -1;  // "EOF"
        }

        // if on new block, seek to its first record
        if (currentBlockReadCount_ == 0) {
            const auto seekResult = bgzf_seek(bgzf, blocks_.at(0).virtualOffset_, SEEK_SET);
            if (seekResult == -1) {
                std::ostringstream s;
                s << "[pbbam] indexed BAM reader  ERROR: could not seek in BAM file:\n"
                  << "  file: " << file_.Filename() << '\n'
                  << "  offset: " << blocks_.at(0).virtualOffset_;
                MaybePrintErrnoReason(s);
                throw std::runtime_error{s.str()};
            }
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

    BamFile file_;
    PbiFilter filter_;
    std::shared_ptr<PbiRawData> index_;
    IndexResultBlocks blocks_;
    size_t currentBlockReadCount_;
    uint32_t numMatchingReads_;
};

PbiIndexedBamReader::PbiIndexedBamReader(PbiFilter filter, const std::string& filename)
    : PbiIndexedBamReader{std::move(filter), BamFile{filename}}
{
}

PbiIndexedBamReader::PbiIndexedBamReader(PbiFilter filter, const std::string& filename,
                                         const std::shared_ptr<PbiRawData>& index)
    : PbiIndexedBamReader{std::move(filter), BamFile{filename}, index}
{
}

PbiIndexedBamReader::PbiIndexedBamReader(PbiFilter filter, BamFile bamFile)
    : PbiIndexedBamReader{std::move(bamFile)}
{
    Filter(std::move(filter));
}

PbiIndexedBamReader::PbiIndexedBamReader(PbiFilter filter, BamFile bamFile,
                                         const std::shared_ptr<PbiRawData>& index)
    : PbiIndexedBamReader{std::move(bamFile), index}
{
    Filter(std::move(filter));
}

PbiIndexedBamReader::PbiIndexedBamReader(const std::string& bamFilename)
    : PbiIndexedBamReader{BamFile{bamFilename}}
{
}

PbiIndexedBamReader::PbiIndexedBamReader(const std::string& bamFilename,
                                         const std::shared_ptr<PbiRawData>& index)
    : PbiIndexedBamReader{BamFile{bamFilename}, index}
{
}

PbiIndexedBamReader::PbiIndexedBamReader(BamFile bamFile) : BamReader{bamFile.Filename()}
{
    auto indexCache = MakePbiIndexCache(bamFile);
    d_ = std::make_unique<PbiIndexedBamReaderPrivate>(std::move(bamFile), indexCache->at(0));
}

PbiIndexedBamReader::PbiIndexedBamReader(BamFile bamFile, const std::shared_ptr<PbiRawData>& index)
    : BamReader{bamFile.Filename()}
    , d_{std::make_unique<PbiIndexedBamReaderPrivate>(std::move(bamFile), index)}
{
}

PbiIndexedBamReader::PbiIndexedBamReader(PbiIndexedBamReader&&) noexcept = default;

PbiIndexedBamReader& PbiIndexedBamReader::operator=(PbiIndexedBamReader&&) noexcept = default;

PbiIndexedBamReader::~PbiIndexedBamReader() = default;

const BamFile& PbiIndexedBamReader::File() const { return d_->file_; }

const PbiFilter& PbiIndexedBamReader::Filter() const { return d_->filter_; }

PbiIndexedBamReader& PbiIndexedBamReader::Filter(PbiFilter filter)
{
    d_->Filter(std::move(filter));
    return *this;
}

const IndexResultBlocks& PbiIndexedBamReader::IndexBlocks() const { return d_->blocks_; }

uint32_t PbiIndexedBamReader::NumReads() const { return d_->numMatchingReads_; }

int PbiIndexedBamReader::ReadRawData(samFile* sf, bam1_t* b)
{
    return d_->ReadRawData(sf->fp.bgzf, b);
}

}  // namespace BAM
}  // namespace PacBio
