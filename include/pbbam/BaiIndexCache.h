#ifndef PBBAM_BAIINDEXCACHE_H
#define PBBAM_BAIINDEXCACHE_H

#include <pbbam/Config.h>

#include <pbcopper/data/Position.h>

#include <htslib/hts.h>

#include <memory>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace BAM {

class BamFile;
class DataSet;

///
/// \brief Caches contents of *.bai file for re-use by multiple readers.
///
class BaiIndexCacheData
{
public:
    explicit BaiIndexCacheData(const BamFile& bamFile);
    explicit BaiIndexCacheData(const std::string& bamFilename);

    BaiIndexCacheData(BaiIndexCacheData&&) noexcept;
    BaiIndexCacheData& operator=(BaiIndexCacheData&&) noexcept;
    ~BaiIndexCacheData();

    /// \note This is very much an internal method and should not be considered
    ///       public API. Exposed here only because of implementation details
    ///       (definition of htslib-related custom deleters) and may be removed.
    ///
    /// \note Does not own the returned pointer; caller is responsible.
    ///
    hts_itr_t* IteratorForInterval(std::int32_t refId, Data::Position start,
                                   Data::Position stop) const;

private:
    struct BaiIndexCacheDataPrivate;
    std::unique_ptr<BaiIndexCacheDataPrivate> d_;
};

using BaiIndexCache = std::shared_ptr<std::vector<std::shared_ptr<BaiIndexCacheData>>>;

BaiIndexCache MakeBaiIndexCache(const DataSet& dataset);
BaiIndexCache MakeBaiIndexCache(const std::vector<BamFile>& bamFiles);
BaiIndexCache MakeBaiIndexCache(const BamFile& bamFile);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BAIINDEXCACHE_H
