// File Description
/// \file BaiIndexedBamReader.cpp
/// \brief Implements the BaiIndexedBamReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BaiIndexedBamReader.h"

#include <cstddef>

#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

struct BaiIndexedBamReaderPrivate
{
public:
    BaiIndexedBamReaderPrivate(const BamFile& file, const GenomicInterval& interval)
        : htsIndex_(nullptr), htsIterator_(nullptr)
    {
        LoadIndex(file.Filename());
        Interval(file.Header(), interval);
    }

    void Interval(const BamHeader& header, const GenomicInterval& interval)
    {
        htsIterator_.reset(nullptr);

        if (header.HasSequence(interval.Name())) {
            auto id = header.SequenceId(interval.Name());
            if (id >= 0 && static_cast<size_t>(id) < header.NumSequences()) {
                htsIterator_.reset(
                    bam_itr_queryi(htsIndex_.get(), id, interval.Start(), interval.Stop()));
            }
        }

        if (!htsIterator_)
            throw std::runtime_error("could not create iterator for requested region");
    }

    void LoadIndex(const std::string& fn)
    {
        htsIndex_.reset(bam_index_load(fn.c_str()));
        if (!htsIndex_) throw std::runtime_error("could not load BAI index data");
    }

    int ReadRawData(BGZF* bgzf, bam1_t* b)
    {
        assert(htsIterator_.get());
        return hts_itr_next(bgzf, htsIterator_.get(), b, nullptr);
    }

public:
    GenomicInterval interval_;
    std::unique_ptr<hts_idx_t, internal::HtslibIndexDeleter> htsIndex_;
    std::unique_ptr<hts_itr_t, internal::HtslibIteratorDeleter> htsIterator_;
};

}  // namespace internal

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval,
                                         const std::string& filename)
    : BaiIndexedBamReader(interval, BamFile(filename))
{
}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, const BamFile& bamFile)
    : BamReader(bamFile), d_(new internal::BaiIndexedBamReaderPrivate(File(), interval))
{
}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, BamFile&& bamFile)
    : BamReader(std::move(bamFile)), d_(new internal::BaiIndexedBamReaderPrivate(File(), interval))
{
}

const GenomicInterval& BaiIndexedBamReader::Interval() const
{
    assert(d_);
    return d_->interval_;
}

int BaiIndexedBamReader::ReadRawData(BGZF* bgzf, bam1_t* b)
{
    assert(d_);
    return d_->ReadRawData(bgzf, b);
}

BaiIndexedBamReader& BaiIndexedBamReader::Interval(const GenomicInterval& interval)
{
    assert(d_);
    d_->Interval(Header(), interval);
    return *this;
}

}  // namespace BAM
}  // namespace PacBio
