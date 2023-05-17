#include "PbbamInternalConfig.h"

#include <pbbam/BaiIndexCache.h>

#include <pbbam/BamFile.h>
#include <pbbam/DataSet.h>
#include <pbbam/Deleters.h>
#include "ErrnoReason.h"

#include <htslib/sam.h>

#include <sstream>
#include <stdexcept>

namespace PacBio {
namespace BAM {

struct BaiIndexCacheData::BaiIndexCacheDataPrivate
{
    using IndexType = std::unique_ptr<hts_idx_t, HtslibIndexDeleter>;
    IndexType htsIndex_;
};

BaiIndexCacheData::BaiIndexCacheData(const BamFile& bamFile) : BaiIndexCacheData(bamFile.Filename())
{}

BaiIndexCacheData::BaiIndexCacheData(const std::string& bamFilename)
    : d_{std::make_unique<BaiIndexCacheData::BaiIndexCacheDataPrivate>()}
{
    d_->htsIndex_.reset(bam_index_load(bamFilename.c_str()));
    if (!d_->htsIndex_) {
        std::ostringstream s;
        s << "[pbbam] BAI index cache ERROR: could not load BAI index data:\n"
          << "  BAM file: " << bamFilename << '\n'
          << "  BAI file: " + bamFilename + ".bai";
        MaybePrintErrnoReason(s);
        throw std::runtime_error{s.str()};
    }
}

BaiIndexCacheData::BaiIndexCacheData(BaiIndexCacheData&&) noexcept = default;

BaiIndexCacheData& BaiIndexCacheData::operator=(BaiIndexCacheData&&) noexcept = default;

BaiIndexCacheData::~BaiIndexCacheData() = default;

hts_itr_t* BaiIndexCacheData::IteratorForInterval(const std::int32_t refId,
                                                  const Data::Position start,
                                                  const Data::Position stop) const
{
    return bam_itr_queryi(d_->htsIndex_.get(), refId, start, stop);
}

using BaiIndexCache = std::shared_ptr<std::vector<std::shared_ptr<BaiIndexCacheData>>>;

BaiIndexCache MakeBaiIndexCache(const DataSet& dataset)
{
    return MakeBaiIndexCache(dataset.BamFiles());
}

BaiIndexCache MakeBaiIndexCache(const std::vector<BamFile>& bamFiles)
{
    auto cache = std::make_shared<std::vector<std::shared_ptr<BaiIndexCacheData>>>();
    auto& indices = *cache.get();
    for (const auto& bamFile : bamFiles) {
        indices.push_back(std::make_shared<BaiIndexCacheData>(bamFile));
    }
    return cache;
}

BaiIndexCache MakeBaiIndexCache(const BamFile& bamFile)
{
    std::vector<BamFile> bamFiles{bamFile};
    return MakeBaiIndexCache(bamFiles);
}

}  // namespace BAM
}  // namespace PacBio
