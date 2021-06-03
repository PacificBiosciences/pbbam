#include "QueryLookup.h"

#include <sstream>
#include <string>

#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/ReadGroupInfo.h>

namespace PacBio {
namespace PbBamify {

QueryLookup::QueryLookup(BAM::DataSet dataset) : dataset_{std::move(dataset)} {}

void QueryLookup::Load()
{
    std::vector<BAM::BamFile> bamFiles{dataset_.BamFiles()};

    // Merge all the read groups for a unified read group lookup.
    BAM::BamHeader jointHeader;
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
        auto new_reader = std::make_shared<BAM::BamReader>(file);
        readers_.push_back(new_reader);
    }

    // Get the PacBio index.
    const BAM::PbiRawData pbi{dataset_};
    const auto& basicData = pbi.BasicData();

    // Clear everything just in case the user called Load() twice.
    lookup_.clear();

    // Process each read in the dataset and reconstruct it's original
    // qname. Place the read in the lookup, together with the ID
    // of the source BAM file and the virtual file offset where
    // the read is located.
    std::ostringstream out;
    for (size_t i = 0; i < pbi.NumReads(); ++i) {
        const auto zmw = basicData.holeNumber_.at(i);
        const auto qStart = basicData.qStart_.at(i);
        const auto qEnd = basicData.qEnd_.at(i);

        const auto rgId = basicData.rgId_.at(i);
        const auto rgString = BAM::ReadGroupInfo::IntToId(rgId);
        const auto rgInfo = jointHeader.ReadGroup(rgString);
        const auto movieName = rgInfo.MovieName();
        std::string type{rgInfo.ReadType()};
        std::transform(type.begin(), type.end(), type.begin(), ::tolower);

        out.str("");
        if (type == "subread") {
            out << movieName << '/' << zmw << '/' << qStart << '_' << qEnd;
        } else if (type == "ccs") {
            out << movieName << '/' << zmw << '/' << "ccs";
        } else {
            out << "Unknown read group type '" << type << "'.";
            throw std::runtime_error(out.str());
        }

        // Sanity check.
        const auto qName = out.str();
        const auto found = lookup_.find(qName);
        if (found != lookup_.end()) {
            const std::string message = std::string{"More than 1 occurrence of qname '"} + qName +
                                        std::string{"'. Duplicate reads in the dataset?"};
            throw std::runtime_error(message);
        }

        const auto fileNumber = basicData.fileNumber_.at(i);
        const auto fileOffset = basicData.fileOffset_.at(i);
        lookup_[qName] = QueryLocation{fileNumber, fileOffset};
    }
}

bool QueryLookup::Find(const std::string& qName, BAM::BamRecord& record) const
{
    const auto it = lookup_.find(qName);
    if (it == lookup_.end()) {
        return false;
    }

    readers_.at(it->second.fileNumber)->VirtualSeek(it->second.fileOffset);
    if (!readers_.at(it->second.fileNumber)->GetNext(record)) {
        return false;
    }

    return true;
}

}  // namespace PbBamify
}  // namespace PacBio
