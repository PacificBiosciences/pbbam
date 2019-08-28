// File Description
/// \file BaiIndexedBamReader.cpp
/// \brief Implements the BaiIndexedBamReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BaiIndexedBamReader.h"

#include <cassert>
#include <cstddef>
#include <sstream>
#include <stdexcept>

#include "MemoryUtils.h"
#include "pbbam/BaiIndexCache.h"

namespace PacBio {
namespace BAM {

class BaiIndexedBamReader::BaiIndexedBamReaderPrivate
{
public:
    BaiIndexedBamReaderPrivate(BamFile file, const std::shared_ptr<BaiIndexCacheData>& index)
        : file_{std::move(file)}, index_{index}
    {
        if (!index_) index_ = std::make_shared<BaiIndexCacheData>(file_);
        assert(index_);  // should throw in cache load if failed
    }

    BaiIndexedBamReaderPrivate(BamFile file, const GenomicInterval& interval,
                               const std::shared_ptr<BaiIndexCacheData>& indexCache)
        : BaiIndexedBamReaderPrivate{std::move(file), indexCache}
    {
        Interval(file_.Header(), interval);
    }

    void Interval(const BamHeader& header, const GenomicInterval& interval)
    {
        htsIterator_.reset();

        if (header.HasSequence(interval.Name())) {
            auto id = header.SequenceId(interval.Name());
            if (id >= 0 && static_cast<size_t>(id) < header.NumSequences()) {
                htsIterator_.reset(
                    index_->IteratorForInterval(id, interval.Start(), interval.Stop()));
            }
        }

        if (!htsIterator_) {
            std::ostringstream s;
            s << "BaiIndexedBamReader: could not create iterator for requested region: "
              << interval.Name() << " [" << interval.Start() << ", " << interval.Stop() << ')';
            throw std::runtime_error{s.str()};
        }
    }

    int ReadRawData(BGZF* bgzf, bam1_t* b)
    {
        assert(htsIterator_.get());
        return hts_itr_next(bgzf, htsIterator_.get(), b, nullptr);
    }

    BamFile file_;
    std::shared_ptr<BaiIndexCacheData> index_;
    GenomicInterval interval_;
    std::unique_ptr<hts_itr_t, HtslibIteratorDeleter> htsIterator_;
};

BaiIndexedBamReader::BaiIndexedBamReader(std::string filename)
    : BaiIndexedBamReader{BamFile{std::move(filename)}, nullptr}
{
}

BaiIndexedBamReader::BaiIndexedBamReader(std::string filename,
                                         const std::shared_ptr<BaiIndexCacheData>& index)
    : BaiIndexedBamReader{BamFile{std::move(filename)}, index}
{
}

BaiIndexedBamReader::BaiIndexedBamReader(BamFile bamFile)
    : BamReader{bamFile.Filename()}
    , d_{std::make_unique<BaiIndexedBamReaderPrivate>(std::move(bamFile), nullptr)}
{
}

BaiIndexedBamReader::BaiIndexedBamReader(BamFile bamFile,
                                         const std::shared_ptr<BaiIndexCacheData>& index)
    : BamReader{bamFile.Filename()}
    , d_{std::make_unique<BaiIndexedBamReaderPrivate>(std::move(bamFile), index)}
{
}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, std::string filename)
    : BaiIndexedBamReader{interval, BamFile{std::move(filename)}, nullptr}
{
}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, std::string filename,
                                         const std::shared_ptr<BaiIndexCacheData>& index)
    : BaiIndexedBamReader{interval, BamFile{std::move(filename)}, index}
{
}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, BamFile bamFile)
    : BamReader{bamFile.Filename()}
    , d_{std::make_unique<BaiIndexedBamReaderPrivate>(std::move(bamFile), interval, nullptr)}
{
}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, BamFile bamFile,
                                         const std::shared_ptr<BaiIndexCacheData>& index)
    : BamReader{bamFile.Filename()}
    , d_{std::make_unique<BaiIndexedBamReaderPrivate>(std::move(bamFile), interval, index)}
{
}

const BamFile& BaiIndexedBamReader::File() const { return d_->file_; }

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
