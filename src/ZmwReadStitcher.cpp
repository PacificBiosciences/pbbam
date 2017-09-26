// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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
/// \file ZmwReadStitcher.cpp
/// \brief Implements the ZmwReadStitcher class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/virtual/ZmwReadStitcher.h"
#include "pbbam/DataSet.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiFilter.h"
#include "pbbam/PbiFilterQuery.h"
#include "VirtualZmwReader.h"
#include <deque>
#include <stdexcept>
#include <utility>

namespace PacBio {
namespace BAM {

struct ZmwReadStitcher::ZmwReadStitcherPrivate
{
public:
    ZmwReadStitcherPrivate(std::string primaryBamFilePath,
                           std::string scrapsBamFilePath,
                           PbiFilter filter)
        : filter_(std::move(filter))
    {
        sources_.push_back(std::make_pair(std::move(primaryBamFilePath), std::move(scrapsBamFilePath)));
        OpenNextReader();
    }

    ZmwReadStitcherPrivate(const DataSet& dataset)
        : filter_(PbiFilter::FromDataSet(dataset))
    {
        // set up source queue
        std::string primaryFn;
        std::string scrapsFn;
        const auto& resources = dataset.ExternalResources();
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
                const auto& childResources = resource.ExternalResources();
                for (const auto& childResource : childResources) {
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
                sources_.emplace_back(make_pair(primaryFn,scrapsFn));
        }

        OpenNextReader();
    }

public:
    bool HasNext() const
    { return (currentReader_ && currentReader_->HasNext()); }

    VirtualZmwBamRecord Next()
    {
        if (currentReader_) {
            const auto result = currentReader_->Next();
            if (!currentReader_->HasNext())
                OpenNextReader();
            return result;
        }

        // no reader active
        const std::string msg = { "no readers active, make sure you use "
                                  "ZmwReadStitcher::HasNext before "
                                  "requesting next record"
                                };
        throw std::runtime_error(msg);
    }

    std::vector<BamRecord> NextRaw()
    {
        if (currentReader_) {
            const auto result = currentReader_->NextRaw();
            if (!currentReader_->HasNext())
                OpenNextReader();
            return result;
        }

        // no reader active
        const std::string msg = { "no readers active, make sure you use "
                                  "ZmwReadStitcher::HasNext before "
                                  "requesting next group of records"
                                };
        throw std::runtime_error(msg);
    }

    BamHeader PrimaryHeader() const
    { return currentReader_->PrimaryHeader(); }

    BamHeader ScrapsHeader() const
    { return currentReader_->ScrapsHeader(); }

private:
    std::deque< std::pair<std::string, std::string> > sources_;
    std::unique_ptr<internal::VirtualZmwReader> currentReader_;
    PbiFilter filter_;

private:
    void OpenNextReader()
    {
        currentReader_.reset(nullptr);

        // find next source pair with data
        while(!sources_.empty()) {
            const auto nextSource = sources_.front();
            sources_.pop_front();

            currentReader_ = std::make_unique<internal::VirtualZmwReader>(nextSource.first,
                                                                          nextSource.second,
                                                                          filter_);
            if (currentReader_->HasNext())
                return;
        }
    }
};

// --------------------------------
// ZmwReadStitcher implementation
// --------------------------------

ZmwReadStitcher::ZmwReadStitcher(const std::string& primaryBamFilePath,
                                 const std::string& scrapsBamFilePath)
    : ZmwReadStitcher(primaryBamFilePath,
                      scrapsBamFilePath,
                      PbiFilter{})
{ }

ZmwReadStitcher::ZmwReadStitcher(const std::string& primaryBamFilePath,
                                 const std::string& scrapsBamFilePath,
                                 const PbiFilter& filter)
    : d_(new ZmwReadStitcherPrivate(primaryBamFilePath,
                                    scrapsBamFilePath,
                                    filter))
{ }

ZmwReadStitcher::ZmwReadStitcher(const DataSet& dataset)
    : d_(new ZmwReadStitcherPrivate(dataset))
{ }

ZmwReadStitcher::~ZmwReadStitcher() { }

bool ZmwReadStitcher::HasNext()
{ return d_->HasNext(); }

VirtualZmwBamRecord ZmwReadStitcher::Next()
{ return d_->Next(); }

std::vector<BamRecord> ZmwReadStitcher::NextRaw()
{ return d_->NextRaw(); }

BamHeader ZmwReadStitcher::PrimaryHeader() const
{ return d_->PrimaryHeader().DeepCopy(); }

BamHeader ZmwReadStitcher::ScrapsHeader() const 
{ return d_->ScrapsHeader().DeepCopy(); }

} // namespace BAM
} // namespace PacBio
