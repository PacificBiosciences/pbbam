#ifndef PBBAM_ZMWCHUNKEDFASTXREADERIMPL_H
#define PBBAM_ZMWCHUNKEDFASTXREADERIMPL_H

#include <pbbam/Config.h>

#include <pbbam/FaiIndex.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/FastqSequence.h>
#include "FaiZmwChunker.h"

#include <string>
#include <utility>

namespace PacBio {
namespace BAM {

class ZmwChunkedFastxReaderImpl
{
public:
    virtual ~ZmwChunkedFastxReaderImpl();

    virtual void Seek(uint64_t pos) = 0;
    virtual FastaSequence ReadNextFasta(bool skipName) = 0;
    virtual FastqSequence ReadNextFastq(bool skipName) = 0;

    std::string fastxFilename_;
    std::string faiFilename_;
    FaiIndex index_;
    FaiZmwChunker chunker_;

protected:
    ZmwChunkedFastxReaderImpl(std::string fastxFilename, std::size_t numChunks);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWCHUNKEDFASTXREADERIMPL_H
