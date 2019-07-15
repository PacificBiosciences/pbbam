// File Description
/// \file IndexedFastqReaderImpl.h
/// \brief Defines the IndexedFastqReaderImpl class.
//
// Author: Derek Barnett

#ifndef INDEXEDFASTQREADERIMPL_H
#define INDEXEDFASTQREADERIMPL_H

#include "pbbam/Config.h"

#include <string>
#include <utility>

#include "pbbam/FaiIndex.h"
#include "pbbam/Position.h"
#include "pbbam/QualityValues.h"

namespace PacBio {
namespace BAM {

class IndexedFastqReaderImpl
{
public:
    virtual ~IndexedFastqReaderImpl();

    virtual std::pair<std::string, QualityValues> Subsequence(const std::string& id, Position start,
                                                              Position end) = 0;

    std::string fastqFilename_;
    std::string faiFilename_;
    FaiIndex index_;

protected:
    IndexedFastqReaderImpl(std::string filename);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // INDEXEDFASTQREADERIMPL_H