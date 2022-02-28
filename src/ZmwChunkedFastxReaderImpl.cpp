#include "PbbamInternalConfig.h"

#include "ZmwChunkedFastxReaderImpl.h"

namespace PacBio {
namespace BAM {

ZmwChunkedFastxReaderImpl::ZmwChunkedFastxReaderImpl(std::string filename, const size_t numChunks)
    : fastxFilename_{std::move(filename)}
    , faiFilename_{fastxFilename_ + ".fai"}
    , index_{faiFilename_}
    , chunker_{index_, numChunks}
{}

ZmwChunkedFastxReaderImpl::~ZmwChunkedFastxReaderImpl() = default;

}  // namespace BAM
}  // namespace PacBio
