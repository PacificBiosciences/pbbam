// Copyright (c) 2017, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Ivan Sovic

#include "QueryLookup.h"

#include <pbbam/PbiRawData.h>
#include <pbbam/PbiLookupData.h>
#include <ostream>
#include <string>
#include <iostream>

namespace PacBio {
namespace BAM {
namespace pbbamify {

std::unique_ptr<QueryLookup> CreateQueryLookup(const PacBio::BAM::DataSet& dataset) {
    return std::unique_ptr<QueryLookup>(new QueryLookup(dataset));
}

QueryLookup::QueryLookup(const PacBio::BAM::DataSet& dataset)
                        : dataset_{dataset} {

}

void QueryLookup::Load() {
    std::vector<BamFile> bamFiles(dataset_.BamFiles());

    // Merge all the read groups for a unified read group lookup.
    PacBio::BAM::BamHeader jointHeader;
    bool headerInitialized = false;
    for (auto& bamFile: bamFiles) {
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
    for (auto& file: bamFiles) {
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
            std::string message = std::string("Unknown read group type '") + type + std::string("'.");
            throw std::runtime_error(message);
        }

        // Sanity check.
        auto it = lookup_.find(qName);
        if (it != lookup_.end()) {
            std::string message = std::string("More than 1 occurrence of qname '") + qName + std::string("'. Duplicate reads in the dataset?");
            throw std::runtime_error(message);
        }

        lookup_[qName] = QueryLocation(fileNumber, fileOffset);
    }
}

bool QueryLookup::Find(const std::string& qName, BamRecord& record) const {
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

} // namespace pbbamify
} // namespace BAM
} // namespace PacBio
