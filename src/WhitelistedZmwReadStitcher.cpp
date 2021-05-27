#include "PbbamInternalConfig.h"

#include <pbbam/virtual/WhitelistedZmwReadStitcher.h>

#include <cstdint>

#include <sstream>
#include <stdexcept>

#include <pbbam/PbiIndexedBamReader.h>

#include "VirtualZmwReader.h"

namespace PacBio {
namespace BAM {

class WhitelistedZmwReadStitcher::WhitelistedZmwReadStitcherPrivate
{
public:
    WhitelistedZmwReadStitcherPrivate(const std::vector<int32_t>& zmwWhitelist,
                                      const std::string& primaryBamFilePath,
                                      const std::string& scrapsBamFilePath)
        : primaryBamFile_{std::make_unique<BamFile>(primaryBamFilePath)}
        , scrapsBamFile_{std::make_unique<BamFile>(scrapsBamFilePath)}
        , primaryReader_{std::make_unique<PbiIndexedBamReader>(*primaryBamFile_)}
        , scrapsReader_{std::make_unique<PbiIndexedBamReader>(*scrapsBamFile_)}
    {
        // setup new header for stitched data
        polyHeader_ = std::make_unique<BamHeader>(primaryBamFile_->Header().ToSam());
        auto readGroups = polyHeader_->ReadGroups();
        if (readGroups.empty()) {
            std::ostringstream msg;
            msg << "[pbbam] ZMW record stitching ERROR: no read groups in header of the primary "
                   "BAM:\n"
                << "  file: " << primaryBamFilePath;
            throw std::runtime_error{msg.str()};
        }
        readGroups[0].ReadType("POLYMERASE");
        readGroups[0].Id(readGroups[0].MovieName(), "POLYMERASE");
        if (readGroups.size() > 1) {
            std::vector<ReadGroupInfo> singleGroup;
            singleGroup.emplace_back(std::move(readGroups[0]));
            readGroups = std::move(singleGroup);
            polyHeader_->ClearReadGroups();
        }
        polyHeader_->ReadGroups(readGroups);

        // remove ZMWs up front, that are not found in either file
        PreFilterZmws(zmwWhitelist);
    }

    bool HasNext() const { return !zmwWhitelist_.empty(); }

    VirtualZmwBamRecord Next()
    {
        auto bamRecordVec = NextRaw();
        return {std::move(bamRecordVec), *polyHeader_};
    }

    std::vector<BamRecord> NextRaw()
    {
        std::vector<BamRecord> result;
        if (!HasNext()) return result;

        const auto& zmw = zmwWhitelist_.front();
        primaryReader_->Filter(PbiZmwFilter{zmw});
        scrapsReader_->Filter(PbiZmwFilter{zmw});

        BamRecord record;
        while (primaryReader_->GetNext(record))
            result.push_back(record);
        while (scrapsReader_->GetNext(record))
            result.push_back(record);

        zmwWhitelist_.pop_front();
        return result;
    }

    BamHeader PrimaryHeader() const { return primaryBamFile_->Header(); }

    BamHeader ScrapsHeader() const { return scrapsBamFile_->Header(); }

private:
    std::unique_ptr<BamFile> primaryBamFile_;
    std::unique_ptr<BamFile> scrapsBamFile_;
    std::unique_ptr<PbiIndexedBamReader> primaryReader_;
    std::unique_ptr<PbiIndexedBamReader> scrapsReader_;
    std::unique_ptr<BamHeader> polyHeader_;
    std::deque<int32_t> zmwWhitelist_;

    void PreFilterZmws(const std::vector<int32_t>& zmwWhitelist)
    {
        // fetch input ZMWs
        const PbiRawData primaryIndex{primaryBamFile_->PacBioIndexFilename()};
        const PbiRawData scrapsIndex{scrapsBamFile_->PacBioIndexFilename()};
        const auto& primaryZmws = primaryIndex.BasicData().holeNumber_;
        const auto& scrapsZmws = scrapsIndex.BasicData().holeNumber_;

        // toss them all into a set (for uniqueness & lookup here soon)
        std::set<int32_t> inputZmws;
        for (const auto& zmw : primaryZmws)
            inputZmws.insert(zmw);
        for (const auto& zmw : scrapsZmws)
            inputZmws.insert(zmw);

        // check our requested whitelist against files' ZMWs, keep if found
        const auto inputEnd = inputZmws.cend();
        for (const int32_t zmw : zmwWhitelist) {
            if (inputZmws.find(zmw) != inputEnd) zmwWhitelist_.push_back(zmw);
        }
    }
};

// --------------------------------
// ZmwReadStitcher implementation
// --------------------------------

WhitelistedZmwReadStitcher::WhitelistedZmwReadStitcher(const std::vector<int32_t>& zmwWhitelist,
                                                       const std::string& primaryBamFilePath,
                                                       const std::string& scrapsBamFilePath)
    : d_{std::make_unique<WhitelistedZmwReadStitcherPrivate>(zmwWhitelist, primaryBamFilePath,
                                                             scrapsBamFilePath)}
{
}

WhitelistedZmwReadStitcher::~WhitelistedZmwReadStitcher() = default;

bool WhitelistedZmwReadStitcher::HasNext() const { return d_->HasNext(); }

VirtualZmwBamRecord WhitelistedZmwReadStitcher::Next() { return d_->Next(); }

std::vector<BamRecord> WhitelistedZmwReadStitcher::NextRaw() { return d_->NextRaw(); }

BamHeader WhitelistedZmwReadStitcher::PrimaryHeader() const { return d_->PrimaryHeader(); }

BamHeader WhitelistedZmwReadStitcher::ScrapsHeader() const { return d_->ScrapsHeader(); }

}  // namespace BAM
}  // namespace PacBio
