// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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
/// \file BamFile.cpp
/// \brief Implements the BamFile class.
//
// Author: Derek Barnett

#include "pbbam/BamFile.h"
#include "pbbam/PbiFile.h"
#include "FileUtils.h"
#include "MemoryUtils.h"
#include <htslib/sam.h>
#include <memory>
#include <sstream>
#include <cassert>
#include <sys/stat.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class BamFilePrivate
{
public:
    BamFilePrivate(const string& fn)
        : filename_(fn)
        , firstAlignmentOffset_(-1)
    {
        // ensure we've updated htslib verbosity with requested verbosity here
        hts_verbose = ( PacBio::BAM::HtslibVerbosity == -1 ? 0 : PacBio::BAM::HtslibVerbosity);

        // attempt open
        auto f = RawOpen();

#if !defined (PBBAM_NO_CHECK_EOF) || PBBAM_AUTOVALIDATE
        // sanity check on file
        const int eofCheck = bgzf_check_EOF(f->fp.bgzf);
        if (eofCheck <= 0 ) {
            // 1:  EOF present & correct
            // 2:  not seekable (e.g. reading from stdin)
            // 0:  EOF absent
            // -1: some other error
            stringstream e;
            if (eofCheck == 0)
                e << fn << " : is missing EOF block" << endl;
            else
                e << fn << " : unknown error while checking EOF block" << endl;
            throw std::runtime_error(e.str());
        }
#endif

        // attempt fetch header
        std::unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> hdr(sam_hdr_read(f.get()));
        header_ = internal::BamHeaderMemory::FromRawData(hdr.get());

        // cache first alignment offset
        firstAlignmentOffset_ = bgzf_tell(f->fp.bgzf);
    }

    unique_ptr<BamFilePrivate> DeepCopy(void)
    {
        return unique_ptr<BamFilePrivate>(new BamFilePrivate(filename_));
    }

    bool HasEOF(void) const
    {
        // streamed input is unknown, since it's not random-accessible
        if (filename_ == "-")
            return false;

        // attempt open
        auto f = RawOpen();
        return RawEOFCheck(f) == 1;
    }

    int RawEOFCheck(const std::unique_ptr<samFile, internal::HtslibFileDeleter>& f) const
    {
        assert(f);
        assert(f->fp.bgzf);
        return bgzf_check_EOF(f->fp.bgzf);
    }

    std::unique_ptr<samFile, internal::HtslibFileDeleter> RawOpen(void) const
    {
        std::unique_ptr<samFile, internal::HtslibFileDeleter> f(sam_open(filename_.c_str(), "rb"));
        if (!f || !f->fp.bgzf)
            throw std::runtime_error(string("could not open BAM file: ") + filename_);
        if (f->format.format != bam)
            throw std::runtime_error("expected BAM, unknown format");
        return f;
    }

public:
    std::string filename_;
    BamHeader header_;
    int64_t firstAlignmentOffset_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

// ------------------------
// BamFile implementation
// ------------------------

BamFile::BamFile(const std::string& filename)
    : d_(new internal::BamFilePrivate(filename))
{ }

BamFile::BamFile(const BamFile& other)
    : d_(other.d_->DeepCopy())
{ }

BamFile::BamFile(BamFile&& other)
    : d_(std::move(other.d_))
{ }

BamFile& BamFile::operator=(const BamFile& other)
{
    d_ = other.d_->DeepCopy();
    return *this;
}

BamFile& BamFile::operator=(BamFile&& other)
{ d_ = std::move(other.d_); return *this; }

BamFile::~BamFile(void) { }

void BamFile::CreatePacBioIndex(void) const
{
    PbiFile::CreateFrom(*this);
}

void BamFile::CreateStandardIndex(void) const
{
    if (bam_index_build(d_->filename_.c_str(), 0) != 0)
        throw std::runtime_error("could not build BAI index");
}

void BamFile::EnsurePacBioIndexExists(void) const
{
    if (!PacBioIndexExists())
        CreatePacBioIndex();
}

void BamFile::EnsureStandardIndexExists(void) const
{
    if (!StandardIndexExists())
        CreateStandardIndex();
}

std::string BamFile::Filename(void) const
{ return d_->filename_; }

int64_t BamFile::FirstAlignmentOffset(void) const
{ return d_->firstAlignmentOffset_; }

bool BamFile::HasEOF(void) const
{ return d_->HasEOF(); }

bool BamFile::HasReference(const std::string& name) const
{ return d_->header_.HasSequence(name); }

const BamHeader& BamFile::Header(void) const
{ return d_->header_; }

bool BamFile::IsPacBioBAM(void) const
{ return !d_->header_.PacBioBamVersion().empty(); }

bool BamFile::PacBioIndexExists(void) const
{ return internal::FileUtils::Exists(PacBioIndexFilename()); }

std::string BamFile::PacBioIndexFilename(void) const
{ return d_->filename_ + ".pbi"; }

bool BamFile::PacBioIndexIsNewer(void) const
{
    const auto bamTimestamp = internal::FileUtils::LastModified(Filename());
    const auto pbiTimestamp = internal::FileUtils::LastModified(PacBioIndexFilename());
    return bamTimestamp <= pbiTimestamp;
}

int BamFile::ReferenceId(const std::string& name) const
{ return d_->header_.SequenceId(name); }

uint32_t BamFile::ReferenceLength(const std::string& name) const
{ return ReferenceLength(ReferenceId(name)); }

uint32_t BamFile::ReferenceLength(const int id) const
{ return std::stoul(d_->header_.SequenceLength(id)); }

std::string BamFile::ReferenceName(const int id) const
{ return d_->header_.SequenceName(id); }

bool BamFile::StandardIndexExists(void) const
{ return internal::FileUtils::Exists(StandardIndexFilename()); }

std::string BamFile::StandardIndexFilename(void) const
{ return d_->filename_ + ".bai"; }

bool BamFile::StandardIndexIsNewer(void) const 
{ 
    const auto bamTimestamp  = internal::FileUtils::LastModified(Filename());
    const auto baiTimestamp = internal::FileUtils::LastModified(StandardIndexFilename());
    return bamTimestamp <= baiTimestamp;
}

