#include "PbbamInternalConfig.h"

#include "VirtualZmwCompositeReader.h"

#include <boost/algorithm/string.hpp>

namespace PacBio {
namespace BAM {

VirtualZmwCompositeReader::VirtualZmwCompositeReader(const DataSet& dataset)
    : currentReader_(nullptr), filter_(PbiFilter::FromDataSet(dataset))
{
    sources_ = SourcesFromDataset(dataset);
    OpenNextReader();
}

bool VirtualZmwCompositeReader::HasNext() { return (currentReader_ && currentReader_->HasNext()); }

VirtualZmwBamRecord VirtualZmwCompositeReader::Next()
{
    if (currentReader_) {
        const auto result = currentReader_->Next();
        if (!currentReader_->HasNext()) OpenNextReader();
        return result;
    }

    // no reader active
    throw std::runtime_error{
        "[pbbam] stitched ZMW record reader ERROR: "
        "no readers active, make sure you use "
        "VirtualZmwCompositeReader::HasNext before "
        "requesting next record"};
}

std::vector<BamRecord> VirtualZmwCompositeReader::NextRaw()
{
    if (currentReader_) {
        const auto result = currentReader_->NextRaw();
        if (!currentReader_->HasNext()) OpenNextReader();
        return result;
    }

    // no reader active
    throw std::runtime_error{
        "[pbbam] stitched ZMW record reader ERROR: "
        "no readers active, make sure you use "
        "VirtualZmwCompositeReader::HasNext before "
        "requesting next group of records"};
}

void VirtualZmwCompositeReader::OpenNextReader()
{
    currentReader_.reset(nullptr);

    // find next source pair with data
    while (!sources_.empty()) {
        const auto nextSource = sources_.front();
        sources_.pop_front();

        currentReader_ =
            std::make_unique<VirtualZmwReader>(nextSource.first, nextSource.second, filter_);
        if (currentReader_->HasNext()) return;
    }
}

}  // namespace BAM
}  // namespace PacBio
