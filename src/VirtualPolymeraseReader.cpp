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
/// \file VirtualPolymeraseReader.cpp
/// \brief Implements the VirtualPolymeraseReader class.
//
// Author: Armin Töpfer

#include <stdexcept>

#include "pbbam/virtual/VirtualPolymeraseReader.h"
#include "pbbam/ReadGroupInfo.h"

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class IBackend
{
protected:
    IBackend(const string& primaryBamFilePath,
             const string& scrapsBamFilePath)
    {
        primaryBamFile_ = std::unique_ptr<BamFile>(new BamFile(primaryBamFilePath));
        scrapsBamFile_  = std::unique_ptr<BamFile>(new BamFile(scrapsBamFilePath));

        polyHeader_     = std::unique_ptr<BamHeader>(
                            new BamHeader(primaryBamFile_->Header().ToSam()));

        auto readGroups = polyHeader_->ReadGroups();
        if (readGroups.empty())
            throw std::runtime_error("Bam header of the primary bam has no read groups.");
        readGroups[0].ReadType("POLYMERASE");
        readGroups[0].Id(readGroups[0].MovieName(), "POLYMERASE");
        if (readGroups.size() > 1)
        {
            std::vector<ReadGroupInfo> singleGroup;
            singleGroup.emplace_back(std::move(readGroups[0]));
            readGroups = std::move(singleGroup);
            polyHeader_->ClearReadGroups();
        }
        polyHeader_->ReadGroups(readGroups);
    }

public:
    ~IBackend(void) { }

public:
    virtual bool HasNext(void) =0;
    virtual std::vector<BamRecord> NextRaw(void) =0;

    const BamHeader& PolyHeader(void) const
    { return *polyHeader_; }

    BamHeader PrimaryHeader(void) const
    { return primaryBamFile_->Header(); }

    BamHeader ScrapsHeader(void) const
    { return scrapsBamFile_->Header(); }

protected:
    std::unique_ptr<BamFile>   primaryBamFile_;
    std::unique_ptr<BamFile>   scrapsBamFile_;
    std::unique_ptr<BamHeader> polyHeader_;
};

class EntireFileBackend : public IBackend
{
public:
    EntireFileBackend(const string& primaryBamFilepath,
                      const string& scrapsBamFilepath)
        : IBackend(primaryBamFilepath, scrapsBamFilepath)
    {
        primaryQuery_   = std::unique_ptr<EntireFileQuery>(new EntireFileQuery(*primaryBamFile_));
        primaryIt_      = primaryQuery_->begin();

        scrapsQuery_    = std::unique_ptr<EntireFileQuery>(new EntireFileQuery(*scrapsBamFile_));
        scrapsIt_       = scrapsQuery_->begin();
    }

    ~EntireFileBackend(void) { }

public:
    bool HasNext(void)
    {
        // Return true until both iterators are at the end of the query
        return primaryIt_ != primaryQuery_->end() || scrapsIt_ != scrapsQuery_->end();
    }

    std::vector<BamRecord> NextRaw(void)
    {
        std::vector<BamRecord> bamRecordVec;

        // Current hole number, the smallest of scraps and primary.
        // It can be that the next ZMW is scrap only.
        int currentHoleNumber;
        if (primaryIt_ == primaryQuery_->end())
            currentHoleNumber = (*scrapsIt_).HoleNumber();
        else if (scrapsIt_ == scrapsQuery_->end())
            currentHoleNumber = (*primaryIt_).HoleNumber();
        else
            currentHoleNumber = std::min((*primaryIt_).HoleNumber(), (*scrapsIt_).HoleNumber());

        // collect subreads or hqregions
        while (primaryIt_ != primaryQuery_->end() && currentHoleNumber == (*primaryIt_).HoleNumber())
            bamRecordVec.push_back(*primaryIt_++);

        // collect scraps
        while (scrapsIt_ != scrapsQuery_->end() && currentHoleNumber == (*scrapsIt_).HoleNumber())
            bamRecordVec.push_back(*scrapsIt_++);

        return bamRecordVec;
    }

    std::unique_ptr<EntireFileQuery> primaryQuery_;
    std::unique_ptr<EntireFileQuery> scrapsQuery_;
    EntireFileQuery::iterator        primaryIt_;
    EntireFileQuery::iterator        scrapsIt_;
};

class PbiFilterBackend : public IBackend
{
public:
    PbiFilterBackend(const string& primaryBamFilePath,
                     const string& scrapsBamFilePath,
                     const PbiFilter& filter)
        : IBackend(primaryBamFilePath, scrapsBamFilePath)
    {
        primaryQuery_   = std::unique_ptr<PbiFilterQuery>(new PbiFilterQuery(filter, *primaryBamFile_));
        primaryIt_      = primaryQuery_->begin();

        scrapsQuery_    = std::unique_ptr<PbiFilterQuery>(new PbiFilterQuery(filter, *scrapsBamFile_));
        scrapsIt_       = scrapsQuery_->begin();
    }

    ~PbiFilterBackend(void) { }

public:
    bool HasNext(void)
    {
        // Return true until both iterators are at the end of the query
        return primaryIt_ != primaryQuery_->end() || scrapsIt_ != scrapsQuery_->end();
    }

    std::vector<BamRecord> NextRaw(void)
    {
        std::vector<BamRecord> bamRecordVec;

        // Current hole number, the smallest of scraps and primary.
        // It can be that the next ZMW is scrap only.
        int currentHoleNumber;
        if (primaryIt_ == primaryQuery_->end())
            currentHoleNumber = (*scrapsIt_).HoleNumber();
        else if (scrapsIt_ == scrapsQuery_->end())
            currentHoleNumber = (*primaryIt_).HoleNumber();
        else
            currentHoleNumber = std::min((*primaryIt_).HoleNumber(), (*scrapsIt_).HoleNumber());

        // collect subreads or hqregions
        while (primaryIt_ != primaryQuery_->end() && currentHoleNumber == (*primaryIt_).HoleNumber())
            bamRecordVec.push_back(*primaryIt_++);

        // collect scraps
        while (scrapsIt_ != scrapsQuery_->end() && currentHoleNumber == (*scrapsIt_).HoleNumber())
            bamRecordVec.push_back(*scrapsIt_++);

        return bamRecordVec;
    }

private:
    std::unique_ptr<PbiFilterQuery> primaryQuery_;
    std::unique_ptr<PbiFilterQuery> scrapsQuery_;
    PbiFilterQuery::iterator        primaryIt_;
    PbiFilterQuery::iterator        scrapsIt_;
};

} // namespace internal

struct VirtualPolymeraseReader::VirtualPolymeraseReaderPrivate
{
    VirtualPolymeraseReaderPrivate(const string& primaryBamFilepath,
                                   const string& scrapsBamFilePath,
                                   const PbiFilter& filter)
        : backend_(nullptr)
    {
        if (filter.IsEmpty()) {
            backend_.reset(new internal::EntireFileBackend(primaryBamFilepath,
                                                           scrapsBamFilePath));
        } else {
            backend_.reset(new internal::PbiFilterBackend(primaryBamFilepath,
                                                          scrapsBamFilePath,
                                                          filter));
        }
    }

    bool HasNext(void)
    { return backend_->HasNext(); }

    std::vector<BamRecord> NextRaw(void)
    { return backend_->NextRaw(); }

    const BamHeader& PolyHeader(void) const
    { return backend_->PolyHeader(); }

    BamHeader PrimaryHeader(void) const
    { return backend_->PrimaryHeader(); }

    BamHeader ScrapsHeader(void) const
    { return backend_->ScrapsHeader(); }

    std::unique_ptr<internal::IBackend> backend_;
};

} // namespace BAM
} // namespace PacBio

VirtualPolymeraseReader::VirtualPolymeraseReader(const std::string& primaryBamFilePath,
                                                 const std::string& scrapsBamFilePath)
    : d_(new VirtualPolymeraseReaderPrivate(primaryBamFilePath, scrapsBamFilePath, PbiFilter()))
{ }

VirtualPolymeraseReader::VirtualPolymeraseReader(const std::string& primaryBamFilePath,
                                                 const std::string& scrapsBamFilePath,
                                                 const PbiFilter& filter)
    : d_(new VirtualPolymeraseReaderPrivate(primaryBamFilePath, scrapsBamFilePath, filter))
{ }

VirtualPolymeraseReader::~VirtualPolymeraseReader(void) { }

bool VirtualPolymeraseReader::HasNext(void)
{ return d_->HasNext(); }

// This method is not thread safe
VirtualPolymeraseBamRecord VirtualPolymeraseReader::Next(void)
{
    auto bamRecordVec = NextRaw();
    VirtualPolymeraseBamRecord stitched(std::move(bamRecordVec), d_->PolyHeader());
    return std::move(stitched);
}

std::vector<BamRecord> VirtualPolymeraseReader::NextRaw(void)
{ return d_->NextRaw(); }

BamHeader VirtualPolymeraseReader::PrimaryHeader(void) const
{ return d_->PrimaryHeader(); }

BamHeader VirtualPolymeraseReader::ScrapsHeader(void) const
{ return d_->ScrapsHeader(); }
