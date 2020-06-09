// File Description
/// \file ZmwChunkedFastxReaderImpl.h
/// \brief Defines the ZmwChunkedFastxReaderImpl class.
//
// Author: Derek Barnett

#ifndef ZMWCHUNKEDFASTXREADERIMPL_H
#define ZMWCHUNKEDFASTXREADERIMPL_H

#include "pbbam/Config.h"

#include <string>
#include <utility>

#include "pbbam/FaiIndex.h"
#include "pbbam/FastaSequence.h"
#include "pbbam/FastqSequence.h"

#include "FaiZmwChunker.h"

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
    ZmwChunkedFastxReaderImpl(std::string fastxFilename, const size_t numChunks);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // ZMWCHUNKEDFASTXREADERIMPL_H
