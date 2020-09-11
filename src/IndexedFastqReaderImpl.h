#ifndef PBBAM_INDEXEDFASTQREADERIMPL_H
#define PBBAM_INDEXEDFASTQREADERIMPL_H

#include <pbbam/Config.h>

#include <string>
#include <utility>

#include <pbcopper/data/Position.h>
#include <pbcopper/data/QualityValues.h>

#include <pbbam/FaiIndex.h>

namespace PacBio {
namespace BAM {

class IndexedFastqReaderImpl
{
public:
    virtual ~IndexedFastqReaderImpl();

    virtual std::pair<std::string, Data::QualityValues> Subsequence(const std::string& id,
                                                                    Data::Position start,
                                                                    Data::Position end) = 0;

    std::string fastqFilename_;
    std::string faiFilename_;
    FaiIndex index_;

protected:
    IndexedFastqReaderImpl(std::string filename);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_INDEXEDFASTQREADERIMPL_H
