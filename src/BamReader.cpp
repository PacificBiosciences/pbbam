// File Description
/// \file BamReader.cpp
/// \brief Implements the BamReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamReader.h"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <sstream>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include "Autovalidate.h"
#include "MemoryUtils.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/Validator.h"

namespace PacBio {
namespace BAM {
namespace internal {

struct BamReaderPrivate
{
public:
    BamReaderPrivate(BamFile bamFile) : bamFile_{std::move(bamFile)} { DoOpen(); }

    void DoOpen()
    {
        // fetch file pointer
        htsFile_.reset(sam_open(bamFile_.Filename().c_str(), "rb"));
        if (!htsFile_) throw std::runtime_error{"could not open BAM file for reading"};
    }

public:
    std::unique_ptr<samFile, internal::HtslibFileDeleter> htsFile_;
    BamFile bamFile_;
};

}  // namespace internal

BamReader::BamReader(std::string fn) : BamReader{BamFile{std::move(fn)}} {}

BamReader::BamReader(BamFile bamFile)
    : d_{std::make_unique<internal::BamReaderPrivate>(std::move(bamFile))}
{
    // skip header
    VirtualSeek(d_->bamFile_.FirstAlignmentOffset());
}

BamReader::~BamReader() {}

BGZF* BamReader::Bgzf() const
{
    assert(d_);
    assert(d_->htsFile_);
    assert(d_->htsFile_->fp.bgzf);
    return d_->htsFile_->fp.bgzf;
}

const BamFile& BamReader::File() const
{
    assert(d_);
    return d_->bamFile_;
}

const std::string& BamReader::Filename() const
{
    assert(d_);
    return d_->bamFile_.Filename();
}

const BamHeader& BamReader::Header() const
{
    assert(d_);
    return d_->bamFile_.Header();
}

bool BamReader::GetNext(BamRecord& record)
{
    assert(Bgzf());
    assert(internal::BamRecordMemory::GetRawData(record).get());

    const auto result = ReadRawData(Bgzf(), internal::BamRecordMemory::GetRawData(record).get());

    // success
    if (result >= 0) {
        internal::BamRecordMemory::UpdateRecordTags(record);
        record.header_ = Header();
        record.ResetCachedPositions();

#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif
        return true;
    }

    // EOF or end-of-data range (not an error)
    else if (result == -1)
        return false;

    // error corrupted file
    else {
        std::ostringstream msg;
        msg << "corrupted BAM file: ";
        if (result == -2)
            msg << "probably truncated";
        else if (result == -3)
            msg << "could not read BAM record's' core data";
        else if (result == -4)
            msg << "could not read BAM record's' variable-length data";
        else
            msg << "unknown reason " + std::to_string(result);
        msg << " (" << Filename() << ')';
        throw std::runtime_error{msg.str()};
    }
}

int BamReader::ReadRawData(BGZF* bgzf, bam1_t* b) { return bam_read1(bgzf, b); }

void BamReader::VirtualSeek(int64_t virtualOffset)
{
    const auto result = bgzf_seek(Bgzf(), virtualOffset, SEEK_SET);
    if (result != 0) throw std::runtime_error{"Failed to seek in BAM file"};
}

int64_t BamReader::VirtualTell() const { return bgzf_tell(Bgzf()); }

}  // namespace BAM
}  // namespace PacBio
