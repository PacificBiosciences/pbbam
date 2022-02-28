#include "PbbamInternalConfig.h"

#include "IndexedFastqReaderImpl.h"

namespace PacBio {
namespace BAM {

IndexedFastqReaderImpl::IndexedFastqReaderImpl(std::string filename)
    : fastqFilename_{std::move(filename)}
    , faiFilename_{fastqFilename_ + ".fai"}
    , index_{faiFilename_}
{}

IndexedFastqReaderImpl::~IndexedFastqReaderImpl() = default;

}  // namespace BAM
}  // namespace PacBio
