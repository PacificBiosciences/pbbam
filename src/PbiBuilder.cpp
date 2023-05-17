#include "PbbamInternalConfig.h"

#include <pbbam/PbiBuilder.h>

#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordImpl.h>
#include <pbbam/Deleters.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/RecordType.h>
#include "ErrnoReason.h"
#include "MemoryUtils.h"
#include "PbiBuilderBase.h"

#include <pbcopper/utility/Deleters.h>

#include <boost/numeric/conversion/cast.hpp>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>

namespace PacBio {
namespace BAM {

class PbiBuilderPrivate : public PacBio::BAM::PbiBuilderBase
{
    // TODO: Make this tweak-able, a la IndexedBamWriter's buffers
    static constexpr std::size_t MAX_BUFFER_SIZE = 0x10000;

public:
    PbiBuilderPrivate(const std::string& pbiFilename, const std::size_t numReferenceSequences,
                      const bool isCoordinateSorted,
                      const PbiBuilder::CompressionLevel compressionLevel,
                      const std::size_t numThreads)
        : PacBio::BAM::PbiBuilderBase{pbiFilename, compressionLevel, numThreads, MAX_BUFFER_SIZE}
    {

        if (isCoordinateSorted && numReferenceSequences > 0) {
            refDataBuilder_ = std::make_unique<PbiReferenceDataBuilder>(numReferenceSequences);
        }
    }
};

// --------------------------------------------
// PbiBuilder - builder API
// --------------------------------------------

PbiBuilder::PbiBuilder(const std::string& pbiFilename, const CompressionLevel compressionLevel,
                       const std::size_t numThreads)
    : PbiBuilder{pbiFilename, 0, false, compressionLevel, numThreads}
{}

PbiBuilder::PbiBuilder(const std::string& pbiFilename, const std::size_t numReferenceSequences,
                       const CompressionLevel compressionLevel, const std::size_t numThreads)
    : PbiBuilder{pbiFilename, numReferenceSequences, (numReferenceSequences > 0), compressionLevel,
                 numThreads}
{}

PbiBuilder::PbiBuilder(const std::string& pbiFilename, const std::size_t numReferenceSequences,
                       const bool isCoordinateSorted, const CompressionLevel compressionLevel,
                       const std::size_t numThreads)
    : d_{std::make_unique<PbiBuilderPrivate>(pbiFilename, numReferenceSequences, isCoordinateSorted,
                                             compressionLevel, numThreads)}
{}

PbiBuilder::PbiBuilder(PbiBuilder&&) noexcept = default;

PbiBuilder& PbiBuilder::operator=(PbiBuilder&&) noexcept = default;

PbiBuilder::~PbiBuilder() noexcept = default;

void PbiBuilder::AddRecord(const BamRecord& record, const int64_t vOffset)
{
    d_->AddRecord(record, vOffset);
}

void PbiBuilder::Close() { d_->Close(); }

}  // namespace BAM
}  // namespace PacBio
