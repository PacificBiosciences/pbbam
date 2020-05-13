// Author: Derek Barnett

#ifndef MEMORYUTILS_H
#define MEMORYUTILS_H

#include "pbbam/Config.h"

#include <memory>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/BamRecordImpl.h"

namespace PacBio {
namespace BAM {

class BamHeaderMemory
{
public:
    static BamHeader FromRawData(bam_hdr_t* header);
    static std::shared_ptr<bam_hdr_t> MakeRawHeader(const BamHeader& header);
};

class BamRecordMemory
{
public:
    static const auto& GetImpl(const BamRecord& r) { return r.impl_; }
    static const auto& GetImpl(const BamRecord* r) { return r->impl_; }

    static const auto& GetRawData(const BamRecordImpl& impl) { return impl.d_; }
    static const auto& GetRawData(const BamRecordImpl* impl) { return impl->d_; }
    static const auto& GetRawData(const BamRecord& r) { return GetRawData(r.impl_); }
    static const auto& GetRawData(const BamRecord* r) { return GetRawData(r->impl_); }

    static void UpdateRecordTags(const BamRecord& r) { UpdateRecordTags(r.impl_); }
    static void UpdateRecordTags(const BamRecordImpl& r) { r.UpdateTagMap(); }
};

}  // namespace BAM
}  // namespace PacBio

#endif  // MEMORYUTILS_H
