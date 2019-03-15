// File Description
/// \file BaiIndexedBamReader.cpp
/// \brief Implements the BaiIndexedBamReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BaiIndexedBamReader.h"

#include <cstddef>
#include <sstream>
#include <stdexcept>

#include "MemoryUtils.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

class BaiIndexedBamReader::BaiIndexedBamReaderPrivate
{
public:
    BaiIndexedBamReaderPrivate(BamFile file) : file_{std::move(file)} { LoadIndex(); }

    BaiIndexedBamReaderPrivate(BamFile file, const GenomicInterval& interval)
        : BaiIndexedBamReaderPrivate{std::move(file)}
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
                    bam_itr_queryi(htsIndex_.get(), id, interval.Start(), interval.Stop()));
            }
        }

        if (!htsIterator_) {
            std::ostringstream s;
            s << "BaiIndexedBamReader: could not create iterator for requested region: "
              << interval.Name() << " [" << interval.Start() << ", " << interval.Stop() << ')';
            throw std::runtime_error{s.str()};
        }
    }

    void LoadIndex()
    {
        const auto& fn = file_.Filename();
        htsIndex_.reset(bam_index_load(fn.c_str()));
        if (!htsIndex_)
            throw std::runtime_error{
                "BaiIndexedBamReader: could not load *.bai index data for file: " + fn};
    }

    int ReadRawData(BGZF* bgzf, bam1_t* b)
    {
        assert(htsIterator_.get());
        return hts_itr_next(bgzf, htsIterator_.get(), b, nullptr);
    }

    BamFile file_;
    GenomicInterval interval_;
    std::unique_ptr<hts_idx_t, HtslibIndexDeleter> htsIndex_;
    std::unique_ptr<hts_itr_t, HtslibIteratorDeleter> htsIterator_;
};

BaiIndexedBamReader::BaiIndexedBamReader(std::string filename)
    : BaiIndexedBamReader{BamFile{std::move(filename)}}
{}

BaiIndexedBamReader::BaiIndexedBamReader(BamFile bamFile)
    : BamReader{bamFile.Filename()}
    , d_{std::make_unique<BaiIndexedBamReaderPrivate>(std::move(bamFile))}
{}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, std::string filename)
    : BaiIndexedBamReader{interval, BamFile{std::move(filename)}}
{}

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval, BamFile bamFile)
    : BamReader{bamFile.Filename()}
    , d_{std::make_unique<BaiIndexedBamReaderPrivate>(std::move(bamFile), interval)}
{}

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
