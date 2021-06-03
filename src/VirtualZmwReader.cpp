#include "PbbamInternalConfig.h"

#include "VirtualZmwReader.h"

#include <sstream>
#include <stdexcept>

#include <pbbam/ReadGroupInfo.h>

namespace PacBio {
namespace BAM {

VirtualZmwReader::VirtualZmwReader(const std::string& primaryBamFilepath,
                                   const std::string& scrapsBamFilepath)
    : VirtualZmwReader(primaryBamFilepath, scrapsBamFilepath, PbiFilter{})
{
}

VirtualZmwReader::VirtualZmwReader(const std::string& primaryBamFilepath,
                                   const std::string& scrapsBamFilepath, const PbiFilter& filter)
{
    primaryBamFile_ = std::make_unique<BamFile>(primaryBamFilepath);
    scrapsBamFile_ = std::make_unique<BamFile>(scrapsBamFilepath);

    if (filter.IsEmpty()) {
        primaryQuery_ = std::make_unique<EntireFileQuery>(*primaryBamFile_);
        scrapsQuery_ = std::make_unique<EntireFileQuery>(*scrapsBamFile_);
    } else {
        primaryQuery_ = std::make_unique<PbiFilterQuery>(filter, *primaryBamFile_);
        scrapsQuery_ = std::make_unique<PbiFilterQuery>(filter, *scrapsBamFile_);
    }

    primaryIt_ = (primaryQuery_->begin());
    scrapsIt_ = (scrapsQuery_->begin());

    stitchedHeader_ = std::make_unique<BamHeader>(primaryBamFile_->Header().ToSam());

    // update stitched read group in header
    auto readGroups = stitchedHeader_->ReadGroups();
    if (readGroups.empty()) {
        std::ostringstream msg;
        msg << "[pbbam] stitched ZMW record reader ERROR: no read groups in header of the primary "
               "BAM:\n"
            << "  file: " << primaryBamFilepath;
        throw std::runtime_error{msg.str()};
    }
    readGroups[0].ReadType("POLYMERASE");
    readGroups[0].Id(readGroups[0].MovieName(), "POLYMERASE");
    if (readGroups.size() > 1) {
        std::vector<ReadGroupInfo> singleGroup;
        singleGroup.emplace_back(std::move(readGroups[0]));
        readGroups = std::move(singleGroup);
        stitchedHeader_->ClearReadGroups();
    }
    stitchedHeader_->ReadGroups(readGroups);
}

VirtualZmwReader::~VirtualZmwReader() = default;

bool VirtualZmwReader::HasNext()
{
    // Return true until both iterators are at the end of the query
    return primaryIt_ != primaryQuery_->end() || scrapsIt_ != scrapsQuery_->end();
}

// This method is not thread safe
VirtualZmwBamRecord VirtualZmwReader::Next()
{
    return VirtualZmwBamRecord{NextRaw(), *stitchedHeader_};
}

std::vector<BamRecord> VirtualZmwReader::NextRaw()
{
    std::vector<BamRecord> bamRecordVec;

    // Current hole number, the smallest of scraps and primary.
    // It can be that the next ZMW is scrap only.
    const int currentHoleNumber = [&]() {
        if (primaryIt_ == primaryQuery_->end()) {
            return (*scrapsIt_).HoleNumber();
        }
        if (scrapsIt_ == scrapsQuery_->end()) {
            return (*primaryIt_).HoleNumber();
        }
        return std::min((*primaryIt_).HoleNumber(), (*scrapsIt_).HoleNumber());
    }();

    // collect subreads or hqregions
    while (primaryIt_ != primaryQuery_->end() && currentHoleNumber == (*primaryIt_).HoleNumber()) {
        bamRecordVec.push_back(*primaryIt_++);
    }

    // collect scraps
    while (scrapsIt_ != scrapsQuery_->end() && currentHoleNumber == (*scrapsIt_).HoleNumber()) {
        bamRecordVec.push_back(*scrapsIt_++);
    }

    return bamRecordVec;
}

BamHeader VirtualZmwReader::PrimaryHeader() const { return primaryBamFile_->Header(); }

BamHeader VirtualZmwReader::ScrapsHeader() const { return scrapsBamFile_->Header(); }

BamHeader VirtualZmwReader::StitchedHeader() const { return *stitchedHeader_.get(); }

}  // namespace BAM
}  // namespace PacBio
