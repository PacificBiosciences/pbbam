// File Description
/// \file BaiIndexCache.cpp
/// \brief Implements the BaiIndexCache class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BaiIndexCache.h"

#include <stdexcept>

#include <htslib/sam.h>

#include "MemoryUtils.h"
#include "pbbam/BamFile.h"
#include "pbbam/DataSet.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

struct BaiIndexCacheData::BaiIndexCacheDataPrivate
{
    using IndexType = std::unique_ptr<hts_idx_t, HtslibIndexDeleter>;
    IndexType htsIndex_;
};

BaiIndexCacheData::BaiIndexCacheData(const BamFile& bamFile) : BaiIndexCacheData(bamFile.Filename())
{
}

BaiIndexCacheData::BaiIndexCacheData(const std::string& bamFilename)
    : d_{std::make_unique<BaiIndexCacheData::BaiIndexCacheDataPrivate>()}
{
    d_->htsIndex_.reset(bam_index_load(bamFilename.c_str()));
    if (!d_->htsIndex_) {
        throw std::runtime_error{"BaiIndexCache: could not load *.bai index data for file: " +
                                 bamFilename};
    }
}

BaiIndexCacheData::~BaiIndexCacheData() = default;

hts_itr_t* BaiIndexCacheData::IteratorForInterval(const int32_t refId, const Position start,
                                                  const Position stop) const
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
    BaiIndexCache cache = std::make_shared<std::vector<std::shared_ptr<BaiIndexCacheData>>>();
    auto& indices = *cache.get();
    for (const auto& bamFile : bamFiles)
        indices.push_back(std::make_shared<BaiIndexCacheData>(bamFile));
    return cache;
}

BaiIndexCache MakeBaiIndexCache(const BamFile& bamFile)
{
    std::vector<BamFile> bamFiles{bamFile};
    return MakeBaiIndexCache(bamFiles);
}

}  // namespace BAM
}  // namespace PacBio
