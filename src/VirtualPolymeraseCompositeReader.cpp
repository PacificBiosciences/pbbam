// Copyright (c) 2015, Pacific Biosciences of California, Inc.
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
//
// File Description
/// \file VirtualPolymeraseCompositeReader.cpp
/// \brief Implements the VirtualPolymeraseCompositeReader class.
//
// Author: Derek Barnett

#include "pbbam/virtual/VirtualPolymeraseCompositeReader.h"
#include <boost/algorithm/string.hpp>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

VirtualPolymeraseCompositeReader::VirtualPolymeraseCompositeReader(const DataSet& dataset)
    : currentReader_(nullptr)
{
    // set up source queue
    string primaryFn;
    string scrapsFn;
    const ExternalResources& resources = dataset.ExternalResources();
    for (const ExternalResource& resource : resources) {

        primaryFn.clear();
        scrapsFn.clear();

        // if resource is possible "primary" BAM
        const auto& metatype = resource.MetaType();
        if (metatype == "PacBio.SubreadFile.SubreadBamFile" ||
            metatype == "PacBio.SubreadFile.HqRegionBamFile")
        {
            // possible resolve relative path
            primaryFn = dataset.ResolvePath(resource.ResourceId());

            // check for associated scraps file
            const ExternalResources& childResources = resource.ExternalResources();
            for (const ExternalResource& childResource : childResources) {
                const auto& childMetatype = childResource.MetaType();
                if (childMetatype == "PacBio.SubreadFile.ScrapsBamFile" ||
                    childMetatype == "PacBio.SubreadFile.HqScrapsBamFile")
                {
                    // possible resolve relative path
                    scrapsFn = dataset.ResolvePath(childResource.ResourceId());
                    break;
                }
            }
        }

        // queue up source for later
        if (!primaryFn.empty() && !scrapsFn.empty())
            sources_.push_back(make_pair(primaryFn,scrapsFn));
    }

    // open first available source
    OpenNextReader();
}

bool VirtualPolymeraseCompositeReader::HasNext(void)
{
    return (currentReader_ && currentReader_->HasNext());
}

VirtualPolymeraseBamRecord VirtualPolymeraseCompositeReader::Next(void)
{
    if (currentReader_) {
        const auto result = currentReader_->Next();
        if (!currentReader_->HasNext())
            OpenNextReader();
        return result;
    }

    // no reader active
    const string msg = { "no readers active, make sure you use "
                         "VirtualPolymeraseCompositeReader::HasNext before "
                         "requesting next record"
                      };
    throw std::runtime_error(msg);
}

vector<BamRecord> VirtualPolymeraseCompositeReader::NextRaw(void)
{
    if (currentReader_) {
        const auto result = currentReader_->NextRaw();
        if (!currentReader_->HasNext())
            OpenNextReader();
        return result;
    }

    // no reader active
    const string msg = { "no readers active, make sure you use "
                         "VirtualPolymeraseCompositeReader::HasNext before "
                         "requesting next group of records"
                      };
    throw std::runtime_error(msg);
}

void VirtualPolymeraseCompositeReader::OpenNextReader(void)
{
    currentReader_.reset(nullptr);

    // find next source pair with data
    while(!sources_.empty()) {
        const auto nextSource = sources_.front();
        sources_.pop_front();

        currentReader_.reset(new VirtualPolymeraseReader(nextSource.first,
                                                         nextSource.second));
        if (currentReader_->HasNext())
            return;
    }
}
