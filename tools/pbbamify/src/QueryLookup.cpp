// Author: Ivan Sovic

#include "QueryLookup.h"

#include <pbbam/PbiRawData.h>
#include <iostream>
#include <ostream>
#include <string>

namespace PacBio {
namespace BAM {
namespace pbbamify {

std::unique_ptr<QueryLookup> CreateQueryLookup(const PacBio::BAM::DataSet& dataset)
{
    return std::unique_ptr<QueryLookup>(new QueryLookup(dataset));
}

QueryLookup::QueryLookup(const PacBio::BAM::DataSet& dataset) : dataset_(dataset) {}

void QueryLookup::Load()
{
    std::vector<BamFile> bamFiles(dataset_.BamFiles());

    // Merge all the read groups for a unified read group lookup.
    PacBio::BAM::BamHeader jointHeader;
    bool headerInitialized = false;
    for (auto& bamFile : bamFiles) {
        auto header = bamFile.Header();
        if (!headerInitialized) {
            jointHeader = header.DeepCopy();
            headerInitialized = true;
        } else {
            jointHeader += header;
        }
    }

    // Set-up a vector of readers for each BAM in the PacBio dataset
    // to allow for random access.
    readers_.clear();
    for (auto& file : bamFiles) {
        auto new_reader = std::make_shared<BamReader>(file);
        readers_.push_back(new_reader);
    }

    // Get the PacBio index.
    PacBio::BAM::PbiRawData pbi(dataset_);
    const auto& basicData = pbi.BasicData();

    // Clear everything just in case the user called Load() twice.
    lookup_.clear();

    // Process each read in the dataset and reconstruct it's original
    // qname. Place the read in the lookup, together with the ID
    // of the source BAM file and the virtual file offset where
    // the read is located.
    for (size_t i = 0; i < pbi.NumReads(); ++i) {
        const auto zmw = basicData.holeNumber_.at(i);
        const auto qStart = basicData.qStart_.at(i);
        const auto qEnd = basicData.qEnd_.at(i);
        const auto& rgId = basicData.rgId_.at(i);
        auto fileNumber = basicData.fileNumber_.at(i);
        auto fileOffset = basicData.fileOffset_.at(i);

        auto rgString = PacBio::BAM::ReadGroupInfo::IntToId(rgId);
        auto rgInfo = jointHeader.ReadGroup(rgString);
        auto type = rgInfo.ReadType();
        std::string movieName = rgInfo.MovieName();

        std::transform(type.begin(), type.end(), type.begin(), ::tolower);

        std::string qName;
        if (type == std::string("subread")) {
            std::ostringstream oss;
            oss << movieName << '/' << zmw << '/' << qStart << '_' << qEnd;
            qName = oss.str();
        } else if (type == std::string("ccs")) {
            std::ostringstream oss;
            oss << movieName << '/' << zmw << '/' << "ccs";
            qName = oss.str();
        } else {
            std::string message =
                std::string("Unknown read group type '") + type + std::string("'.");
            throw std::runtime_error(message);
        }

        // Sanity check.
        auto it = lookup_.find(qName);
        if (it != lookup_.end()) {
            std::string message = std::string("More than 1 occurrence of qname '") + qName +
                                  std::string("'. Duplicate reads in the dataset?");
            throw std::runtime_error(message);
        }

        lookup_[qName] = QueryLocation(fileNumber, fileOffset);
    }
}

bool QueryLookup::Find(const std::string& qName, BamRecord& record) const
{
    auto it = lookup_.find(qName);

    if (it == lookup_.end()) {
        return false;
    }

    readers_.at(it->second.fileNumber)->VirtualSeek(it->second.fileOffset);

    if (!readers_.at(it->second.fileNumber)->GetNext(record)) {
        return false;
    }

    return true;
}

}  // namespace pbbamify
}  // namespace BAM
}  // namespace PacBio
